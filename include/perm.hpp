/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_HPP
#define PERM_HPP

#include "perm_common_types.hpp"
#include <cstddef> // For size_t
#include <functional>
#include <ostream>

namespace perm {
struct parameters_in_t {
    size_t steps = 1000;
    size_t monomers = 100;
    /// Used only if num_monomers = 0
    float_t end_to_end_distance = 0.0;
    void print(std::ostream &os) const;
};

struct parameters_out_t {
    parameters_out_t() = default;
    parameters_out_t(const parameters_in_t &input) : in(input){};
    single_chain_t<int> chain;
    parameters_in_t in;
    float_t energy = 0.0;
    void print(std::ostream &os) const;
};

single_chain_t<int> random_walk_lattice_2D(const size_t &monomers);
parameters_out_t run_simple_sampling(const parameters_in_t &parameters);
float_t
energy(const single_chain_t<int> &chain,
       const std::function<float_t(const vec3D_t<int> &, const vec3D_t<int> &)>
               &energy_pair_func);

} // namespace perm
#endif
