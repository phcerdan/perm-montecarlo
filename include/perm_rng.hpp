/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_RNG_HPP
#define PERM_RNG_HPP

#include "lattice_lut.hpp"
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
    std::uniform_int_distribution<int> uid(min, max);
    return uid(engine());
}

inline perm::vec3D_t<int> rand_lattice_1D_2n() {
    static thread_local std::uniform_int_distribution<int> uid(0, 1);
    return perm::lattice_1D_2n.at(uid(engine()));
}

inline perm::vec3D_t<int> rand_lattice_2D_4n() {
    static thread_local std::uniform_int_distribution<int> uid(0, 3);
    return perm::lattice_2D_4n.at(uid(engine()));
}

inline perm::vec3D_t<int> rand_lattice_2D_8n() {
    static thread_local std::uniform_int_distribution<int> uid(0, 7);
    return perm::lattice_2D_8n.at(uid(engine()));
}

inline perm::vec3D_t<int> rand_lattice_3D_6n() {
    static thread_local std::uniform_int_distribution<int> uid(0, 5);
    return perm::lattice_3D_6n.at(uid(engine()));
}
inline perm::vec3D_t<int> rand_lattice_3D_18n() {
    static thread_local std::uniform_int_distribution<int> uid(0, 17);
    return perm::lattice_3D_18n.at(uid(engine()));
}
inline perm::vec3D_t<int> rand_lattice_3D_26n() {
    static thread_local std::uniform_int_distribution<int> uid(0, 25);
    return perm::lattice_3D_26n.at(uid(engine()));
}

} // namespace RNG

#endif
