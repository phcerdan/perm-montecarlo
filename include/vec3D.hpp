/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_VEC3D_HPP
#define PERM_VEC3D_HPP

#include "perm_common_types.hpp"
#include <cmath>
#include <iostream>

namespace perm {

template <typename T>
struct vec3D_t {
    T x,y,z;
    void print(std::ostream &os) const {
        os << x << " " << y << " " << z << std::endl;
    }
    bool operator==(const vec3D_t<T> &other) const {
        return x == other.x && y == other.y && z == other.z;
    }
    T &operator[](const size_t index) {
        return (&x)[index];
    }
    const T &operator[](const size_t index) const {
        return (&x)[index];
    }
    constexpr size_t size() { return 3; };
};

struct vec3D_int_hasher {
    std::size_t operator()(const vec3D_t<int> &vec) const {
        auto vec_copy = vec;
        std::size_t h = 0;
        h ^= std::hash<int>{}(vec_copy.x)  + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(vec_copy.y)  + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(vec_copy.z)  + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
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
} // end namespace perm
    
#endif
