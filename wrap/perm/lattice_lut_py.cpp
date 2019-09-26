/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "lattice_lut.hpp"
#include "perm_common_py.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace perm;

void init_lattice_lut(py::module &m) {
    py::bind_map<perm::lattice_map_t>(m, "lattice_map_t");
    py::module mlattice = m.def_submodule(
            "lattice", "Functions dealing with the lattice look-up tables");
    mlattice.attr("d2_n4") = py::cast(lattice_2D_4n);
    mlattice.attr("d2_n8") = py::cast(lattice_2D_8n);
    mlattice.attr("d3_n6") = py::cast(lattice_3D_6n);
    mlattice.attr("d3_n18") = py::cast(lattice_3D_18n);
    mlattice.attr("d3_n26") = py::cast(lattice_3D_26n);
    // py::class_<dumb_lattice>(m, "lattice")
    //         .def_property_readonly_static(
    //                 "d2_n4",
    //                 [](py::object /* self */) { return lattice_2D_4n; })
    //         .def_property_readonly_static(
    //                 "d2_n8",
    //                 [](py::object /* self */) { return lattice_2D_8n; })
    //         .def_property_readonly_static(
    //                 "d3_n6",
    //                 [](py::object /* self */) { return lattice_3D_6n; })
    //         .def_property_readonly_static(
    //                 "d3_n18",
    //                 [](py::object /* self */) { return lattice_3D_18n; })
    //         .def_property_readonly_static("d3_n26", [](py::object /* self */) {
    //             return lattice_3D_26n;
    //         });
}
