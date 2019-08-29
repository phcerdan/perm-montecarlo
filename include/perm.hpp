/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_HPP
#define PERM_HPP

#include "perm_common_types.hpp"
#include <cmath>
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

single_chain_t<int>
random_walk_lattice(const size_t &monomers,
                    const std::function<vec3D_t<int>(void)> &rand_lattice_func);
single_chain_t<int> random_walk_lattice(const size_t &monomers,
                                        const size_t &dimension,
                                        const size_t &neighbors);
parameters_out_t run_simple_sampling(const parameters_in_t &parameters);
float_t
energy(const single_chain_t<int> &chain,
       const std::function<float_t(const vec3D_t<int> &, const vec3D_t<int> &)>
               &energy_pair_func);

template <typename T>
vec3D_t<float_t> end_to_end_vector(const single_chain_t<T> &chain) {
    if (chain.points.empty()) {
        return vec3D_t<float_t>();
    }
    const auto start = chain.points[0];
    const auto end = chain.points.back();
    return {static_cast<float_t>(end.x - start.x),
            static_cast<float_t>(end.y - start.y),
            static_cast<float_t>(end.z - start.z)};
}

template <typename T>
float_t end_to_end_distance(const single_chain_t<T> &chain) {
    if (chain.points.empty()) {
        return 0;
    }
    return distance(chain.points[0], chain.points.back());
}
} // namespace perm
#endif
