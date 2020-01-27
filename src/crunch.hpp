#pragma once

// xtensor includes
#include "pybind11/pybind11.h"
#include "xtensor/xmath.hpp"
#include "xtensor/xarray.hpp"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pyvectorize.hpp"
// end xtensor includes

// xt::pyarray<double> resolution_independent_random_grid(std::size_t gridsize, unsigned int seed);
// xt::pyarray<double> naive_random_grid(std::size_t gridsize, unsigned int seed);
