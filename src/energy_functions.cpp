/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "energy_functions.hpp"
#include <cassert>

namespace perm {

perm::float_t chain_total_energy(const single_chain_t<int> &input_chain,
                                 energy_grow_func_t &energy_func) {
    single_chain_t<int> growing_chain;
    perm::float_t total_energy = 0.0;
    for (const auto &monomer : input_chain.points) {
        total_energy += energy_func(growing_chain, monomer);
        growing_chain.points.push_back(monomer);
        growing_chain.monomers++;
    }
    assert(growing_chain.points.size() == input_chain.points.size());
    return total_energy;
}

perm::float_t energy_grow_zero(const single_chain_t<int> & /*input_chain*/,
                               const vec3D_t<int> & /*new_monomer*/) {
    return 0.0;
}
perm::float_t energy_grow_bending(const single_chain_t<int> &input_chain,
                                  const vec3D_t<int> &new_monomer,
                                  const perm::float_t &k) {
    const auto chain_size = input_chain.points.size();
    // This energy needs at least 2 beads in the chain
    if (chain_size <= 1) {
        return 0.0;
    }
    const auto &a = input_chain.points[chain_size - 2];
    const auto &middle = input_chain.points.back();
    const auto &c = new_monomer;
    return k * (1.0 - perm::cos_director(perm::minus(middle, a),
                                         perm::minus(c, middle)));
}
} // namespace perm
