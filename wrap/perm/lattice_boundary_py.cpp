/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "lattice_boundary.hpp"
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace perm;

void init_lattice_boundary(py::module &m) {
    py::module mboundary = m.def_submodule(
            "boundary", "Functions dealing with volume boundaries");
    mboundary.def("is_inside_always", &boundary::is_inside_always);
    mboundary.def(
            "is_inside_hyper_cube",
            [](const size_t &half_side) -> boundary_func_t {
                return [half_side = half_side](const vec3D_t<int> &point) {
                    return boundary::is_inside_hyper_cube<3>(point, half_side);
                };
            },
            py::arg("half_side"));
    mboundary.def(
            "is_inside_hyper_rectangle",
            [](const size_t &half_side_x, const size_t &half_side_y,
               const size_t &half_side_z) -> boundary_func_t {
                return [half_side_x = half_side_x, half_side_y = half_side_y,
                        half_side_z = half_side_z](const vec3D_t<int> &point) {
                    return boundary::is_inside_hyper_rectangle<3>(
                            point, half_side_x, half_side_y, half_side_z);
                };
            },
            py::arg("half_side_x"), py::arg("half_side_y"),
            py::arg("half_side_z"));
}
