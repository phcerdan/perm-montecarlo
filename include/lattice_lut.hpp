/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_LATTICE_LUT_HPP
#define PERM_LATTICE_LUT_HPP

#include "perm_common_types.hpp"
#include "vec3D.hpp"
#include <unordered_map>
namespace perm {

/// Correct average Kuhn monomer length for second nearest neighbors lattices.
inline float_t bond_length_lattice(const size_t &dim, const size_t &neighbors) {
    float_t bmean = 1.0;

    if (neighbors == 8) {
        bmean = (4 * 1 + 4 * sqrt(2)) / neighbors;
    } else if (neighbors == 18) {
        bmean = (6 * 1 + 12 * sqrt(2)) / neighbors;
    } else if (neighbors == 26) {
        bmean = (6 * 1 + 12 * sqrt(2) + 8 * sqrt(3)) / neighbors;
    }
    return bmean;
}

// clang-format off
using lattice_map_t = std::unordered_map<int, perm::vec3D_t<int>>;
/// Maps are sorted
//

const std::unordered_map<int, perm::vec3D_t<int>> lattice_1D_2n = {
        {0, perm::vec3D_t<int>{-1, 0, 0}},
        {1, perm::vec3D_t<int>{ 1, 0, 0}}};
const std::unordered_map<int, int> inverse_direction_lattice_1D_2n = {
        {0, 1},
        {1, 0}};

const std::unordered_map<int, perm::vec3D_t<int>> lattice_2D_4n = {
        {0, perm::vec3D_t<int>{-1,  0, 0}},
        {1, perm::vec3D_t<int>{ 0, -1, 0}},
        {2, perm::vec3D_t<int>{ 0,  1, 0}},
        {3, perm::vec3D_t<int>{ 1,  0, 0}}};
const std::unordered_map<int, int> inverse_direction_lattice_2D_4n = {
        {0, 3},
        {1, 2},
        {2, 1},
        {3, 0}};

const std::unordered_map<int, perm::vec3D_t<int>> lattice_2D_8n = {
        {0, perm::vec3D_t<int>{-1, -1, 0}},
        {1, perm::vec3D_t<int>{-1,  0, 0}},
        {2, perm::vec3D_t<int>{-1,  1, 0}},
        {3, perm::vec3D_t<int>{ 0, -1, 0}},
        {4, perm::vec3D_t<int>{ 0,  1, 0}},
        {5, perm::vec3D_t<int>{ 1, -1, 0}},
        {6, perm::vec3D_t<int>{ 1,  0, 0}},
        {7, perm::vec3D_t<int>{ 1,  1, 0}}};
const std::unordered_map<int, int> inverse_direction_lattice_2D_8n = {
        {0, 7},
        {1, 6},
        {2, 5},
        {3, 4},
        {4, 3},
        {5, 2},
        {6, 1},
        {7, 0}};

const std::unordered_map<int, perm::vec3D_t<int>> lattice_3D_6n = {
        {0, perm::vec3D_t<int>{-1,  0,  0}},
        {1, perm::vec3D_t<int>{ 0, -1,  0}},
        {2, perm::vec3D_t<int>{ 0,  0, -1}},
        {3, perm::vec3D_t<int>{ 0,  0,  1}},
        {4, perm::vec3D_t<int>{ 0,  1,  0}},
        {5, perm::vec3D_t<int>{ 1,  0,  0}}};
const std::unordered_map<int, int> inverse_direction_lattice_3D_6n = {
        {0, 5},
        {1, 4},
        {2, 3},
        {3, 2},
        {4, 1},
        {5, 0}};

const std::unordered_map<int, perm::vec3D_t<int>> lattice_3D_18n = {
        {0 , perm::vec3D_t<int>{-1, -1,  0}},
        {1 , perm::vec3D_t<int>{-1,  0, -1}},
        {2 , perm::vec3D_t<int>{-1,  0,  0}},
        {3 , perm::vec3D_t<int>{-1,  0,  1}},
        {4 , perm::vec3D_t<int>{-1,  1,  0}},
        {5 , perm::vec3D_t<int>{ 0, -1, -1}},
        {6 , perm::vec3D_t<int>{ 0, -1,  0}},
        {7 , perm::vec3D_t<int>{ 0, -1,  1}},
        {8 , perm::vec3D_t<int>{ 0,  0, -1}},
        {9 , perm::vec3D_t<int>{ 0,  0,  1}},
        {10, perm::vec3D_t<int>{ 0,  1, -1}},
        {11, perm::vec3D_t<int>{ 0,  1,  0}},
        {12, perm::vec3D_t<int>{ 0,  1,  1}},
        {13, perm::vec3D_t<int>{ 1, -1,  0}},
        {14, perm::vec3D_t<int>{ 1,  0, -1}},
        {15, perm::vec3D_t<int>{ 1,  0,  0}},
        {16, perm::vec3D_t<int>{ 1,  0,  1}},
        {17, perm::vec3D_t<int>{ 1,  1,  0}}};
const std::unordered_map<int, int> inverse_direction_lattice_3D_18n = {
        {0 , 17},
        {1 , 16},
        {2 , 15},
        {3 , 14},
        {4 , 13},
        {5 , 12},
        {6 , 11},
        {7 , 10},
        {8 ,  9},
        {9 ,  8},
        {10,  7},
        {11,  6},
        {12,  5},
        {13,  4},
        {14,  3},
        {15,  2},
        {16,  1},
        {17,  0}};

const std::unordered_map<int, perm::vec3D_t<int>> lattice_3D_26n = {
        {0 , perm::vec3D_t<int>{-1, -1, -1}},
        {1 , perm::vec3D_t<int>{-1, -1,  0}},
        {2 , perm::vec3D_t<int>{-1, -1,  1}},
        {3 , perm::vec3D_t<int>{-1,  0, -1}},
        {4 , perm::vec3D_t<int>{-1,  0,  0}},
        {5 , perm::vec3D_t<int>{-1,  0,  1}},
        {6 , perm::vec3D_t<int>{-1,  1, -1}},
        {7 , perm::vec3D_t<int>{-1,  1,  0}},
        {8 , perm::vec3D_t<int>{-1,  1,  1}},
        {9 , perm::vec3D_t<int>{ 0, -1, -1}},
        {10, perm::vec3D_t<int>{ 0, -1,  0}},
        {11, perm::vec3D_t<int>{ 0, -1,  1}},
        {12, perm::vec3D_t<int>{ 0,  0, -1}},
        {13, perm::vec3D_t<int>{ 0,  0,  1}},
        {14, perm::vec3D_t<int>{ 0,  1, -1}},
        {15, perm::vec3D_t<int>{ 0,  1,  0}},
        {16, perm::vec3D_t<int>{ 0,  1,  1}},
        {17, perm::vec3D_t<int>{ 1, -1, -1}},
        {18, perm::vec3D_t<int>{ 1, -1,  0}},
        {19, perm::vec3D_t<int>{ 1, -1,  1}},
        {20, perm::vec3D_t<int>{ 1,  0, -1}},
        {21, perm::vec3D_t<int>{ 1,  0,  0}},
        {22, perm::vec3D_t<int>{ 1,  0,  1}},
        {23, perm::vec3D_t<int>{ 1,  1, -1}},
        {24, perm::vec3D_t<int>{ 1,  1,  0}},
        {25, perm::vec3D_t<int>{ 1,  1,  1}}};

const std::unordered_map<int, int> inverse_direction_lattice_3D_26n = {
        {0 , 25},
        {1 , 24},
        {2 , 23},
        {3 , 22},
        {4 , 21},
        {5 , 20},
        {6 , 19},
        {7 , 18},
        {8 , 17},
        {9 , 16},
        {10, 15},
        {11, 14},
        {12, 13},
        {13, 12},
        {14, 11},
        {15, 10},
        {16,  9},
        {17,  8},
        {18,  7},
        {19,  6},
        {20,  5},
        {21,  4},
        {22,  3},
        {23,  2},
        {24,  1},
        {25,  0}};

// clang-format on
inline const std::unordered_map<int, int>
get_inverse_direction_lattice(const lattice_map_t &lattice) {
    if (lattice == lattice_1D_2n) {
        return inverse_direction_lattice_1D_2n;
    } else if (lattice == lattice_2D_4n) {
        return inverse_direction_lattice_2D_4n;
    } else if (lattice == lattice_2D_8n) {
        return inverse_direction_lattice_2D_8n;
    } else if (lattice == lattice_3D_6n) {
        return inverse_direction_lattice_3D_6n;
    } else if (lattice == lattice_3D_18n) {
        return inverse_direction_lattice_3D_18n;
    } else if (lattice == lattice_3D_26n) {
        return inverse_direction_lattice_3D_26n;
    } else {
        throw std::runtime_error(
                "inverse_direction_lattice for lattice not known");
    }
}
} // namespace perm

#endif
