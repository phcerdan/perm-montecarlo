/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "lattice_lut.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace perm;

void init_lattice_lut(py::module &m) {
    py::bind_map<std::unordered_map<int, vec3D_t<int>>>(m, "lattice_map_t");
    py::class_<perm::lattice_map_t>(m, "lattice")
            .def_property_readonly_static(
                    "d2_n4",
                    [](py::object /* self */) { return lattice_2D_4n; })
            .def_property_readonly_static(
                    "d2_n8",
                    [](py::object /* self */) { return lattice_2D_8n; })
            .def_property_readonly_static(
                    "d3_n6",
                    [](py::object /* self */) { return lattice_3D_6n; })
            .def_property_readonly_static(
                    "d3_n18",
                    [](py::object /* self */) { return lattice_3D_18n; })
            .def_property_readonly_static("d3_n26", [](py::object /* self */) {
                return lattice_3D_26n;
            });
}
