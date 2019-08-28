/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "perm.hpp"
#include "perm_rng.hpp"
#include <iostream>

namespace perm {

void parameters_in_t::print(std::ostream &os) const {
    os << "steps= " << steps << std::endl;
    os << "monomers= " << monomers << std::endl;
    os << "end_to_end_distance= " << end_to_end_distance << std::endl;
}

void parameters_out_t::print(std::ostream &os) const {
    in.print(os);
    os << "energy= " << energy << std::endl;
    chain.print(os);
}

parameters_out_t run_simple_sampling(const parameters_in_t &parameters_in) {
    parameters_out_t parameters_out(parameters_in);
    auto &chain = parameters_out.chain;
    chain.points.emplace_back(vec3D_t<int>{0, 0, 0});
    chain.monomers++;
    return parameters_out;
}

single_chain_t<int> random_walk_lattice_2D(const size_t &monomers) {
    single_chain_t<int> chain;
    chain.points.emplace_back(vec3D_t<int>{0, 0, 0});
    chain.monomers++;
    while (chain.monomers < monomers) {
        chain.points.emplace_back(
                perm::plus(chain.points.back(), RNG::rand_lattice_2D()));
        chain.monomers++;
    }
    return chain;
}

float_t
energy(const single_chain_t<int> &chain,
       const std::function<float_t(const vec3D_t<int> &, const vec3D_t<int> &)>
               &energy_pair_func) {

    if (chain.monomers == 0 || chain.points.empty()) {
        return 0.0;
    }

    auto &points = chain.points;
    double sum = 0.0;
    for (size_t i = 0; i != points.size() - 1; i++) {
        sum += energy_pair_func(points[i], points[i + 1]);
    }
    return sum;
}
} // namespace perm
