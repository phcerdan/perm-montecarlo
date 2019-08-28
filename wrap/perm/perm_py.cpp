/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "perm.hpp"
#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace perm;

void init_perm(py::module &m) {
    m.def("random_walk_lattice_2D", &random_walk_lattice_2D);
}
