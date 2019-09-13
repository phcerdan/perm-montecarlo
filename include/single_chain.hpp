/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_SINGLE_CHAIN_HPP
#define PERM_SINGLE_CHAIN_HPP

#include "vec3D.hpp"
#include <vector>

namespace perm {

template <typename T>
struct single_chain_t {
    /// Number of monomers in the single-chain polymer
    size_t monomers = 0;
    /** Ordered collection of points,
     * start: points[0],
     * end: points[num_monomers - 1]
     */
    std::vector<vec3D_t<T>> points;
    void print(std::ostream &os) const {
        os << "chain.monomers= " << monomers << std::endl;
        for (const auto &p : points) {
            p.print(os);
        }
    }
};

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

template <typename T>
float_t contour_length(const single_chain_t<T> &chain) {
    const auto num_monomers = chain.points.size();
    if (num_monomers == 0 || num_monomers == 1) {
        return 0;
    }

    float_t dist = 0.0;
    for (size_t i = 0; i < num_monomers - 1; i++) {
        dist += perm::distance(chain.points[i], chain.points[i + 1]);
    }
    return dist;
}

template <typename T>
vec3D_t<float_t> center_of_mass(const single_chain_t<T> &chain) {
    if (chain.points.empty()) {
        return vec3D_t<float_t>();
    }
    auto vec_center = vec3D_t<float_t>();
    const auto num_monomers = chain.points.size();
    for (int i = 0; i < num_monomers; i++) {
        vec_center = perm::plus<float_t, T, float_t>(
                vec_center, perm::minus(chain.points[i], chain.points[0]));
    }
    return perm::multiply(vec_center, 1.0 / num_monomers);
}

template <typename T>
float_t gyration_radius_square(const single_chain_t<T> &chain) {
    const auto center = center_of_mass(chain);
    float_t out = 0.0;
    const auto num_monomers = chain.points.size();
    for (int i = 0; i < num_monomers; i++) {
        out += perm::distance_square(center, chain.points[i]);
    }
    return out / num_monomers;
}
} // namespace perm
#endif
