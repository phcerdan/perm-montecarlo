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

single_chain_t<int> random_walk_lattice(
        const size_t &monomers,
        const std::function<vec3D_t<int>(void)> &rand_lattice_func) {
    single_chain_t<int> chain;
    chain.points.emplace_back(vec3D_t<int>{0, 0, 0});
    chain.monomers++;
    while (chain.monomers < monomers) {
        chain.points.emplace_back(
                perm::plus(chain.points.back(), rand_lattice_func()));
        chain.monomers++;
    }
    return chain;
}
single_chain_t<int> random_walk_lattice(const size_t &monomers,
                                        const size_t &dimension,
                                        const size_t &neighbors) {
    if (dimension == 1) {
        if (neighbors == 2) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_1D_2n);
        } else {
            throw std::logic_error(
                    "In 1D, neighbors only be 2, but neighbors= " +
                    std::to_string(neighbors));
        }
    } else if (dimension == 2) {
        if (neighbors == 4) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_2D_4n);
        } else if (neighbors == 8) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_2D_8n);
        } else {
            throw std::logic_error(
                    "In 2D, neighbors can be 4 or 8 but neighbors= " +
                    std::to_string(neighbors));
        }
    } else if (dimension == 3) {
        if (neighbors == 6) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_3D_6n);
        } else if (neighbors == 18) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_3D_18n);
        } else if (neighbors == 26) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_3D_26n);
        } else {
            throw std::logic_error(
                    "In 3D, neighbors can be 6, 18 or 26 but neighbors= " +
                    std::to_string(neighbors));
        }
    } else {
        throw std::logic_error(
                "Only dimension 1, 2 or 3 supported but dimension= " +
                std::to_string(dimension));
    }
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
