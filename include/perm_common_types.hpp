/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_COMMON_TYPES_HPP
#define PERM_COMMON_TYPES_HPP
#include <iostream>
#include <vector>

namespace perm {
/// float type (hard-coded better than templated for now)
using float_t = double;

template <typename T>
struct vec3D_t {
    T x;
    T y;
    T z;
    void print(std::ostream &os) const {
        os << x << " " << y << " " << z << std::endl;
    }
    friend bool operator==(const vec3D_t<T> &lhs, const vec3D_t<T> &rhs) {
        return lhs.x == lhs.x && lhs.y == lhs.y && lhs.z == lhs.z;
    }
    T &operator[](const size_t index) {
        return (index == 0 ? x : (index == 1 ? y : z));
    }
    const T &operator[](const size_t index) const {
        return (index == 0 ? x : (index == 1 ? y : z));
    }
    constexpr size_t size() { return 3; };
};

template <typename T>
vec3D_t<T> plus(const vec3D_t<T> &lhs, const vec3D_t<T> &rhs) {
    vec3D_t<T> out;
    out.x = lhs.x + rhs.x;
    out.y = lhs.y + rhs.y;
    out.z = lhs.z + rhs.z;
    return out;
}

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
} // namespace perm
#endif
