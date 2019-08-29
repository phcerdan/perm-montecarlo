/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "declare_types_py.hpp"
#include "perm.hpp"
#include "perm_common_py.hpp"
#include "perm_common_types.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>

namespace py = pybind11;
using namespace perm;

void init_perm_common_types(py::module &m) {
    perm::declare_vec3D<int>(m, "i");
    perm::declare_vec3D<perm::float_t>(m, "f");
    py::bind_vector<std::vector<vec3D_t<int>>>(m, "points_vector",
                                               py::module_local(false));
    py::class_<single_chain_t<int>>(m, "single_chain")
            .def(py::init())
            .def_readwrite("monomers", &single_chain_t<int>::monomers)
            .def_readwrite("points", &single_chain_t<int>::points)
            .def("ete_vector",
                 [](const single_chain_t<int> &chain) {
                     return end_to_end_vector(chain);
                 })
            .def("ete_distance",
                 [](const single_chain_t<int> &chain) {
                     return end_to_end_distance(chain);
                 })
            .def("center_of_mass",
                 [](const single_chain_t<int> &chain) {
                     return center_of_mass(chain);
                 })
            .def("gyration_radius_square",
                 [](const single_chain_t<int> &chain) {
                     return gyration_radius_square(chain);
                 })
            .def("__repr__", [](const single_chain_t<int> &chain) {
                std::stringstream os;
                chain.print(os);
                return os.str();
            });
}
