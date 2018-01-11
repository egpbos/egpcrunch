#include "pybind11/pybind11.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xarray.hpp"

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include <iostream>
#include <numeric>
#include <cmath>
#include <exception>
#include <string>

namespace py = pybind11;

// density estimators

using ULONG = unsigned long;
using real_prec = double;
using ndarray = xt::pyarray<real_prec>;

real_prec pow_3(real_prec a) {
  return a * a * a;
}

real_prec SPH_kernel_3D(real_prec r, real_prec h) {
  // N.B.: when using this for density estimation, you need to normalize the
  // result afterwards with V/N! See e.g. density_SPH. Same goes for
  // grad_SPH_kernel_3D.

  // Monaghan kernel W_4
  real_prec result = 0.;
  real_prec q = r/h;
  if (q < 0.) {
    std::string message = "in SPH_kernel_3D: q is negative!";
    throw std::runtime_error(message.c_str());
  } else if (q <= 1.) {
    result = 1./M_PI/(h*h*h) * (1 - 3./2*q*q + 3./4*q*q*q);
  } else if (q <= 2.) {
    result = 1./M_PI/(h*h*h) * (1./4 * pow_3(2.-q));
  }

  return(result);
}

ndarray density_SPH(ULONG N1, ULONG N2, ULONG N3,
                    real_prec L1, real_prec L2, real_prec L3,
                    real_prec d1, real_prec d2, real_prec d3,
                    real_prec min1, real_prec min2, real_prec min3,
                    const ndarray &xp, const ndarray &yp,
                    const ndarray &zp, const ndarray &masses,
                    bool use_masses, real_prec kernel_h) {
  // Notes: - this implementation assumes periodic boundary conditions!
  //.       - this function returns rho, i.e. not a normalized overdensity
  //
  // Warning: don't use with kernel_h larger than N{1,2,3}/4, otherwise the
  // kernel diameter will be larger than the box width in that direction!
  // When using this function inside Barcode using the supplied INIT_PARAMS
  // initialization function, this condition will always be enforced.
  // 
  // Warning 2: Make sure N{1,2,3} are smaller than LONG_MAX, since they need to
  // be cast to signed integers.

  xt::pyarray<real_prec> delta({N1, N2, N3}, 0);

  ULONG N_OBJ = xp.size();
  if (xp.size() != yp.size() || yp.size() != zp.size() || (use_masses && zp.size() != masses.size())) {
    std::string message = "in density_SPH: particle array lengths do not match!";
    throw std::runtime_error(message.c_str());
  }

  // First determine the reach of the kernel, which determines how many cells
  // every particle needs to loop over to check for contributions. This is
  // based on the kernel radius ( == 2*kernel_h ).
  int reach1 = static_cast<int>(2*kernel_h/d1) + 1;
  int reach2 = static_cast<int>(2*kernel_h/d2) + 1;
  int reach3 = static_cast<int>(2*kernel_h/d3) + 1;

  ULONG NLOSS=0;
  #pragma omp parallel for
  for (ULONG n = 0; n < N_OBJ; ++n) {
    // check if particle is in selected domain, else discard it
    if((xp[n] >= min1 && xp[n] < min1+L1) &&
       (yp[n] >= min2 && yp[n] < min2+L2) &&
       (zp[n] >= min3 && zp[n] < min3+L3)) {
      real_prec mass = 1;
      if (use_masses) {
        mass = masses[n];
      }

      // Determine central cell index where particle resides
      ULONG ix = static_cast<ULONG>(xp[n]/d1);
      ULONG iy = static_cast<ULONG>(yp[n]/d2);
      ULONG iz = static_cast<ULONG>(zp[n]/d3);
      // Central cell position:
      real_prec ccx = (static_cast<real_prec>(ix) + 0.5)*d1;
      real_prec ccy = (static_cast<real_prec>(iy) + 0.5)*d2;
      real_prec ccz = (static_cast<real_prec>(iz) + 0.5)*d3;

      // Loop over surrounding gridcells (including central cell itself) within
      // kernel radius.
      for(int i1 = -reach1; i1 <= reach1; ++i1) {
        for(int i2 = -reach2; i2 <= reach2; ++i2) {
          for(int i3 = -reach3; i3 <= reach3; ++i3) {
            // Cell position (relative to the central cell):
            real_prec cx = ccx + static_cast<real_prec>(i1)*d1;
            real_prec cy = ccy + static_cast<real_prec>(i2)*d2;
            real_prec cz = ccz + static_cast<real_prec>(i3)*d3;
            // Cell index, taking into account periodic boundary conditions:
            ULONG kx, ky, kz;
            kx = (static_cast<ULONG>(static_cast<long>(N1) + i1) + ix) % N1;  // checked in testcodes/signed_unsigned_periodic.cpp that signedness casts here doesn't cause bugs
            ky = (static_cast<ULONG>(static_cast<long>(N2) + i2) + iy) % N2;
            kz = (static_cast<ULONG>(static_cast<long>(N3) + i3) + iz) % N3;
            // The above casts are necessary to ensure implicit signedness casts don't cause trouble.
            // We assume here that kernel_h (and thus reach) <= N{1,2,3}/4 and N1 < LONG_MAX (the latter is true in Barcode, since Ni are usually uint, not ulong).

            real_prec diff_x = xp[n] - cx;
            real_prec diff_y = yp[n] - cy;
            real_prec diff_z = zp[n] - cz;
            real_prec r = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
            
            if (r/kernel_h <= 2.) {
              #pragma omp atomic
              delta(kx, ky, kz) += SPH_kernel_3D(r, kernel_h) * mass;
              // TODO: check whether this is the right order of indices!
            }
          }
        }
      }
    } else {
      ++NLOSS;
    }
  }
  if (NLOSS>0) {
    std::cout << "Lost " << NLOSS << " out of " << N_OBJ << " particles in density_SPH." << std::endl;
  }

  return delta;
}


// version for simple density boxes:
// cubic box, cubic grid cells, no range stuff and no masses
ndarray density_SPH_simple(ULONG Nx, real_prec Lx,
                           const ndarray &xp, const ndarray &yp,
                           const ndarray &zp, real_prec kernel_h) {
  ndarray masses;
  real_prec dx = Lx / static_cast<real_prec>(Nx);
  return density_SPH(Nx, Nx, Nx, Lx, Lx, Lx, dx, dx, dx, 0, 0, 0,
                     xp, yp, zp, masses, false, kernel_h);
}


// Python Module and Docstrings

PYBIND11_MODULE(egpcrunch, m) {
    xt::import_numpy();

    m.doc() = R"pbdoc(
        Number crunching module, used in cosmology research.

        .. currentmodule:: egpcrunch

        .. autosummary::
           :toctree: _generate

           density_SPH
           density_SPH_simple
    )pbdoc";

    m.def("density_SPH", density_SPH, "Estimate the density of a set of particles using fixed-size SPH kernels");
    m.def("density_SPH_simple", density_SPH_simple, "Simple version of density_SPH, cubic box, no range selection and no masses");
}
