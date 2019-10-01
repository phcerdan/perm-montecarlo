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
#include <pybind11/functional.h> // for funcs in  parameters_in_t
#include <sstream>

namespace py = pybind11;
using namespace perm;

void init_perm_common_types(py::module &m) {
    perm::declare_vec3D<int>(m, "i");
    perm::declare_vec3D<perm::float_t>(m, "f");
    py::bind_vector<std::vector<vec3D_t<int>>>(m, "vector_points",
                                               py::module_local(false));
    py::bind_vector<std::vector<perm::float_t>>(m, "vector_float");
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
            .def("contour_length",
                 [](const single_chain_t<int> &chain) {
                     return contour_length(chain);
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
    py::class_<parameters_in_t>(m, "parameters_in_t")
            .def(py::init())
            .def(py::init([](const size_t &num_monomers) {
                return parameters_in_t(num_monomers);
            }))
            .def_readwrite("max_tries", &parameters_in_t::max_tries)
            .def_readwrite("monomers", &parameters_in_t::monomers)
            .def_readwrite("weight_threshold_high",
                           &parameters_in_t::weight_threshold_high)
            .def_readwrite("weight_threshold_low",
                           &parameters_in_t::weight_threshold_low)
            .def_readwrite("lattice", &parameters_in_t::lattice)
            .def_readwrite("end_to_end_distance",
                           &parameters_in_t::end_to_end_distance)
            .def_readwrite("energy_grow_func",
                           &parameters_in_t::energy_grow_func)
            .def_readwrite("is_inside_boundary_func",
                           &parameters_in_t::is_inside_boundary_func)
            .def("__repr__", [](const parameters_in_t &param) {
                std::stringstream os;
                param.print(os);
                return os.str();
                });
    py::class_<parameters_out_many_chains_t>(m, "parameters_out_many_chains_t")
            .def(py::init())
            .def_readwrite("num_chains", &parameters_out_many_chains_t::num_chains)
            .def_readwrite("chains", &parameters_out_many_chains_t::chains)
            .def_readwrite("weights", &parameters_out_many_chains_t::weights)
            .def("__repr__", [](const parameters_out_many_chains_t &param) {
                std::stringstream os;
                param.print(os);
                return os.str();
                });
}
