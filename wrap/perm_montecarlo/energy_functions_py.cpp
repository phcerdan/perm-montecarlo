/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "energy_functions.hpp"
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace perm;

void init_energy_functions(py::module &m) {
    m.def("chain_total_energy", &chain_total_energy);
    m.def("energy_grow_zero", &energy_grow_zero);
    m.def(
            "energy_grow_bending",
            [](const perm::float_t &k) -> energy_grow_func_t {
                return [k = k](const single_chain_t<int> &input_chain,
                               const vec3D_t<int> &new_monomer) {
                    return energy_grow_bending(input_chain, new_monomer, k);
                };
            },
            py::arg("k"));
}
