#include <iostream>
#include <numeric>
#include <cmath>
#include <exception>
#include <string>
#include <random>

#include "crunch.hpp"

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



// random grid creation

xt::pyarray<double> resolution_independent_random_grid(std::size_t gridsize, unsigned int seed) {
    std::mt19937 engine; // only initialize an RNG once! http://tinyurl.com/cwbqqvg
    
    // Create the output array:
    xt::pyarray<double> out({gridsize, gridsize, (gridsize / 2 + 1)}, 0);

    // Seed random number generator:
    engine.seed(seed); // don't seed too much! http://tinyurl.com/cwbqqvg
    
    // Fill out:
    for (std::size_t i = 0; i < gridsize / 2; i++) {
        for (std::size_t k = 0; k < i+1; k++) {
            for (std::size_t j = 0; j < i; j++)
                // outpoint[i * i_stride + j * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out(i, j, k) = (double) ((unsigned long) engine()) / 0x100000000;
            for (std::size_t j = 0; j < i + 1; j++)
                // outpoint[j * i_stride + i * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out(j, i, k) = (double) ((unsigned long) engine()) / 0x100000000;
            for (std::size_t j = 0; j < i; j++)
                // outpoint[(gridsize - 1 - i) * i_stride + j * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out((gridsize - 1 - i), j, k) = (double) ((unsigned long) engine()) / 0x100000000;
            for (std::size_t j = 0; j < i + 1; j++)
                // outpoint[(gridsize - 1 - j) * i_stride + i * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out((gridsize - 1 - j), i, k) = (double) ((unsigned long) engine()) / 0x100000000;
            for (std::size_t j = 0; j < i; j++)
                // outpoint[i * i_stride + (gridsize - 1 - j) * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out(i, (gridsize - 1 - j), k) = (double) ((unsigned long) engine()) / 0x100000000;
            for (std::size_t j = 0; j < i + 1; j++)
                // outpoint[j * i_stride + (gridsize - 1 - i) * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out(j, (gridsize - 1 - i), k) = (double) ((unsigned long) engine()) / 0x100000000;
            for (std::size_t j = 0; j < i; j++)
                // outpoint[(gridsize - 1 - i) * i_stride + (gridsize - 1 - j) * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out((gridsize - 1 - i), (gridsize - 1 - j), k) = (double) ((unsigned long) engine()) / 0x100000000;
            for (std::size_t j = 0; j < i + 1; j++)
                // outpoint[(gridsize - 1 - j) * i_stride + (gridsize - 1 - i) * j_stride + k] = (double) ((unsigned long) engine()) / 0x100000000;
                out((gridsize - 1 - j), (gridsize - 1 - i), k) = (double) ((unsigned long) engine()) / 0x100000000;
        }
        for (std::size_t j = 0; j < i; j++) {
            for (std::size_t k = 0; k < i; k++) {
                // outpoint[j * i_stride + k * j_stride + i] = (double) ((unsigned long) engine()) / 0x100000000;
                out(j, k, i) = (double) ((unsigned long) engine()) / 0x100000000;
                // outpoint[(gridsize - 1 - j) * i_stride + k * j_stride + i] = (double) ((unsigned long) engine()) / 0x100000000;
                out((gridsize - 1 - j), k, i) = (double) ((unsigned long) engine()) / 0x100000000;
                // outpoint[j * i_stride + (gridsize - 1 - k) * j_stride + i] = (double) ((unsigned long) engine()) / 0x100000000;
                out(j, (gridsize - 1 - k), i) = (double) ((unsigned long) engine()) / 0x100000000;
                // outpoint[(gridsize - 1 - j) * i_stride + (gridsize - 1 - k) * j_stride + i] = (double) ((unsigned long) engine()) / 0x100000000;
                out((gridsize - 1 - j), (gridsize - 1 - k), i) = (double) ((unsigned long) engine()) / 0x100000000;
            }
        }
    }
    
    // Fill out's nyquist plane:
    for (std::size_t i = 0; i < gridsize; i++) {
        for (std::size_t j = 0; j < gridsize; j++) {
            // outpoint[i_stride*i + j_stride*j + gridsize/2] = (double) ((unsigned long) engine()) / 0x100000000;
            out(i, j, gridsize/2) = (double) ((unsigned long) engine()) / 0x100000000;
        }
    }
    
    return out;
}

// EXTRA TEST FUNCTION FOR COMPARISON WITH numpy.random.random
xt::pyarray<double> naive_random_grid(std::size_t gridsize, unsigned int seed) {
    std::mt19937 engine; // only initialize an RNG once! http://tinyurl.com/cwbqqvg
    
    // Create the output array:
    xt::pyarray<double> out({gridsize, gridsize, (gridsize/2+1)}, 0);

    // Seed random number generator:
    engine.seed(seed); // don't seed too much! http://tinyurl.com/cwbqqvg
    
    // Fill out:
    for (std::size_t i = 0; i < gridsize / 2; i++) {
        for (std::size_t k = 0; k < gridsize; k++) {
            for (std::size_t j = 0; j < gridsize; j++) {
                out(i, j, k) = (double) ((unsigned long) engine()) / 0x100000000;
            }
        }
    }
    
    return out;
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
           resolution_independent_random_grid
           naive_random_grid
    )pbdoc";

    m.def("density_SPH", density_SPH, "Estimate the density of a set of particles using fixed-size SPH kernels");
    m.def("density_SPH_simple", density_SPH_simple, "Simple version of density_SPH, cubic box, no range selection and no masses");
    m.def("resolution_independent_random_grid", resolution_independent_random_grid, "Build a random grid that can be used to build resolution independent random fields");
    m.def("naive_random_grid", naive_random_grid, "Build a random grid that can NOT be used to build resolution independent random fields, i.e. the same way as regular numpy.random.random");
}
