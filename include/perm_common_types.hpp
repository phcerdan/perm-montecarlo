/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_COMMON_TYPES_HPP
#define PERM_COMMON_TYPES_HPP
#include <cmath>
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

template <typename TLhs, typename TRhs = TLhs, typename TOut = TLhs>
vec3D_t<TOut> plus(const vec3D_t<TLhs> &lhs, const vec3D_t<TRhs> &rhs) {
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

template <typename TLhs, typename TRhs = TLhs, typename TOut = TLhs>
vec3D_t<TOut> minus(const vec3D_t<TLhs> &lhs, const vec3D_t<TRhs> &rhs) {
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

template <typename TLhs, typename TScalar = TLhs, typename TOut = TLhs>
vec3D_t<TOut> multiply(const vec3D_t<TLhs> &lhs, const TScalar &scalar) {
    return {lhs.x * scalar, lhs.y * scalar, lhs.z * scalar};
}

template <typename T>
vec3D_t<T> negate(const vec3D_t<T> &lhs) {
    return {-lhs.x, -lhs.y, -lhs.z};
}

template <typename TLhs, typename TRhs = TLhs>
perm::float_t dot_product(const vec3D_t<TLhs> &lhs, const vec3D_t<TRhs> &rhs) {
    return static_cast<perm::float_t>(lhs.x * rhs.x + lhs.y * rhs.y +
                                      lhs.z * rhs.z);
}
template <typename TLhs, typename TRhs = TLhs, typename TOut = TLhs>
vec3D_t<TOut> cross_product(const vec3D_t<TLhs> &lhs,
                            const vec3D_t<TRhs> &rhs) {
    return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
}

template <typename T>
perm::float_t norm(const vec3D_t<T> &lhs) {
    return sqrt(dot_product(lhs, lhs));
}

template <typename TLhs, typename TRhs = TLhs>
perm::float_t distance_square(const vec3D_t<TLhs> &lhs,
                              const vec3D_t<TRhs> &rhs) {
    const auto diff = minus(lhs, rhs);
    return dot_product(diff, diff);
}

template <typename TLhs, typename TRhs = TLhs>
perm::float_t distance(const vec3D_t<TLhs> &lhs, const vec3D_t<TRhs> &rhs) {
    return norm(minus(lhs, rhs));
}

template <typename TLhs, typename TRhs = TLhs>
perm::float_t angle(const vec3D_t<TLhs> &lhs, const vec3D_t<TRhs> &rhs) {
    return std::atan2(norm(cross_product(lhs, rhs)), dot_product(lhs, rhs));
}

template <typename TLhs, typename TRhs = TLhs>
perm::float_t cos_director(const vec3D_t<TLhs> &lhs, const vec3D_t<TRhs> &rhs) {
    return std::cos(angle(lhs, rhs));
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
