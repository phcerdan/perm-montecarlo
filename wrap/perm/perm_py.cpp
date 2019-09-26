/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "energy_functions.hpp"
#include "lattice_lut.hpp"
#include "perm.hpp"
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace perm;

void init_perm(py::module &m) {
    m.def(
            "random_walk_lattice",
            [](const size_t &monomers, const size_t &dimension,
               const size_t &neighbors

            ) { return random_walk_lattice(monomers, dimension, neighbors); },
            py::arg("N"), py::arg("dim"), py::arg("neighbors"));
    m.def("bond_length_lattice", &perm::bond_length_lattice, py::arg("dim"),
          py::arg("neighbors"));
    m.def("mc_saw_simple_sampling", &perm::mc_saw_simple_sampling,
          py::arg("monomers"), py::arg("tries"), py::arg("lattice"));
    m.def("mc_saw_rosenbluth_sampling", &perm::mc_saw_rosenbluth_sampling,
          py::arg("monomers"), py::arg("tries"), py::arg("lattice"));
    m.def(
            "mc_saw_perm",
            [](const size_t &monomers, const size_t &tries,
               const lattice_map_t &lattice,
               const energy_grow_func_t &energy_grow_func,
               const boundary_func_t &is_inside_boundary_func,
               const std::vector<perm::float_t> &weight_threshold_low,
               const std::vector<perm::float_t> &weight_threshold_high) {
                parameters_in_t parameters_in(monomers);
                parameters_in.max_tries = tries;
                parameters_in.lattice = lattice;
                parameters_in.energy_grow_func = energy_grow_func;
                parameters_in.is_inside_boundary_func = is_inside_boundary_func;
                if (!weight_threshold_low.empty()) {
                    if (weight_threshold_low.size() != monomers) {
                        throw std::runtime_error(
                                "vector weight_threshold_low must have "
                                "number of components equal to monomers in the "
                                "chain");
                    }
                    parameters_in.weight_threshold_low = weight_threshold_low;
                }
                if (!weight_threshold_high.empty()) {
                    if (weight_threshold_high.size() != monomers) {
                        throw std::runtime_error(
                                "vector weight_threshold_high must have "
                                "number of components equal to monomers in the "
                                "chain");
                    }
                    parameters_in.weight_threshold_high = weight_threshold_high;
                }
                return perm::mc_saw_perm(parameters_in);
            },
            py::arg("monomers"), py::arg("tries"), py::arg("lattice"),
            py::arg("energy_grow_func"), py::arg("is_inside_boundary_func"),
            // [](const vec3D_t<int> &) -> bool {return true;},
            py::arg("weight_threshold_low") = std::vector<perm::float_t>(),
            py::arg("weight_threshold_high") = std::vector<perm::float_t>());
    m.def(
            "generate_chains",
            [](const size_t &num_chains, const parameters_in_t &parameters_in,
               const std::string &method) {
        if (method == "perm") {
            return generate_chains(num_chains, mc_saw_perm, parameters_in);
            } else if (method == "rosenbluth") {
                return generate_chains(num_chains, mc_saw_rosenbluth_sampling,
                                       parameters_in.monomers,
                                       parameters_in.max_tries,
                                       parameters_in.lattice);
            }
            else {
                throw std::runtime_error("invalid method, try 'perm'");
            }
        },
            py::arg("num_chains"), py::arg("parameters_in"),
            py::arg("method") = "perm");
}
