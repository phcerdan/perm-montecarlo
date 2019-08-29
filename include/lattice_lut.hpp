/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_LATTICE_LUT_HPP
#define PERM_LATTICE_LUT_HPP

#include "perm_common_types.hpp"
#include <unordered_map>
namespace perm {

// clang-format off
/// Maps are sorted
const std::unordered_map<int, perm::vec3D_t<int>> lattice_1D_2n = {
        {0, perm::vec3D_t<int>{-1, 0, 0}},
        {1, perm::vec3D_t<int>{ 1, 0, 0}},
};
const std::unordered_map<int, perm::vec3D_t<int>> lattice_2D_4n = {
        {0, perm::vec3D_t<int>{-1,  0, 0}},
        {1, perm::vec3D_t<int>{ 0, -1, 0}},
        {2, perm::vec3D_t<int>{ 0,  1, 0}},
        {3, perm::vec3D_t<int>{ 1,  0, 0}}};

const std::unordered_map<int, perm::vec3D_t<int>> lattice_2D_8n = {
        {0, perm::vec3D_t<int>{-1, -1, 0}},
        {1, perm::vec3D_t<int>{-1,  0, 0}},
        {2, perm::vec3D_t<int>{-1,  1, 0}},
        {3, perm::vec3D_t<int>{ 0, -1, 0}},
        {4, perm::vec3D_t<int>{ 0,  1, 0}},
        {5, perm::vec3D_t<int>{ 1, -1, 0}},
        {6, perm::vec3D_t<int>{ 1,  0, 0}},
        {7, perm::vec3D_t<int>{ 1,  1, 0}}};

const std::unordered_map<int, perm::vec3D_t<int>> lattice_3D_6n = {
        {0, perm::vec3D_t<int>{-1,  0,  0}},
        {1, perm::vec3D_t<int>{ 0, -1,  0}},
        {2, perm::vec3D_t<int>{ 0,  0, -1}},
        {3, perm::vec3D_t<int>{ 0,  0,  1}},
        {4, perm::vec3D_t<int>{ 0,  1,  0}},
        {5, perm::vec3D_t<int>{ 1,  0,  0}}};

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
} // namespace perm
// clang-format on

#endif
