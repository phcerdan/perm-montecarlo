/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "perm_common_types.hpp"
#include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
// for PYBIND11_NUMPY_DTYPE(perm::vec3D_t<int>, x, y, z);
// to pass numpy arrays to c++. UNUSED for now

namespace py = pybind11;
void init_perm_common_types(py::module &);
void init_perm(py::module &);
void init_lattice_lut(py::module &);
void init_energy_functions(py::module &);

// using vec3D_int_t = perm::vec3D_t<int>;

PYBIND11_MODULE(_perm, m) {
    m.doc() = "_perm"; // optional module docstring
    init_lattice_lut(m);
    init_perm_common_types(m);
    init_energy_functions(m);
    init_perm(m);
}
