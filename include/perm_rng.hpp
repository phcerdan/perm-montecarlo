/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_RNG_HPP
#define PERM_RNG_HPP

#include "perm_common_types.hpp"
#include <random>
namespace RNG {

inline std::mt19937 &engine() {
    /// seed generation
    static thread_local std::random_device rdev{};
    /// engine instantiation: e
    static thread_local std::mt19937 e{rdev()};
    // static thread_local std::mt19937 e{4342};
    return e;
}

/// Randomize reseed the input engine with a random generate seed.
inline void randomize_engine(std::mt19937 &eng) {
    std::random_device rd{};
    eng.seed(rd());
}

/// Uniform random distribution from double [0.0,1)
inline double rand01() {
    static thread_local std::uniform_real_distribution<double> uid(0.0, 1.0);
    return uid(engine());
}

/**
 * Return 1 with probability p
 * @param p must be lesser or equal than 1
 * @return 1 with probability p, 0 if not.
 */
inline bool random_bool(const double p) { return (rand01() < p) ? 1 : 0; }

/**
 *  Uniform random distribution from int [min,max]
 * @param min
 * @param max
 * @return int from [min,max]
 */
inline int rand_range_int(const int &min, const int &max) {
    // note that inside function static variables doesn't interfer if they have
    // the same name
    std::uniform_int_distribution<int> uid(min, max);
    return uid(engine());
}

inline perm::vec3D_t<int> rand_lattice_2D() {
    // note that inside function static variables doesn't interfer if they have
    // the same name
    std::uniform_int_distribution<int> uid(0, 3);
    const auto lattice_int = uid(engine());
    switch (lattice_int) {
    case 0 /* -x */:
        return perm::vec3D_t<int>{-1, 0, 0};
        break;
    case 1 /* +x */:
        return perm::vec3D_t<int>{1, 0, 0};
        break;
    case 2 /* -y */:
        return perm::vec3D_t<int>{0, -1, 0};
        break;
    case 3 /* +y */:
        return perm::vec3D_t<int>{0, 1, 0};
        break;
    }
    // Not really needed, unreachable, but warnings if removed in gcc
    return perm::vec3D_t<int>{999, 999, 999};
}
} // namespace RNG

#endif
