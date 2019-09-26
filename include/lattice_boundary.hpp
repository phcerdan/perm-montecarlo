/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_LATTICE_BOUNDARY_HPP
#define PERM_LATTICE_BOUNDARY_HPP
#include "vec3D.hpp"
#include <cmath> // for abs
#include <functional>

namespace perm {

using boundary_func_t = std::function<bool(const vec3D_t<int> &)>;

namespace boundary {
constexpr bool is_inside = true;
constexpr bool is_outside = false;
inline bool is_inside_always(const vec3D_t<int> &point) { return is_inside; };

/**
 * Check if point is inside rectangle.
 * Assume hyper_rectangle is centered in 0,0,0
 *
 * @tparam TDimension dimension of the hyper rectangle. options: 1, 2 or 3
 * @param point input point to check if in or out boundaries
 * @param parameter_half_length_x half side of the rectangle in x
 * @param parameter_half_length_y half side of the rectangle in y
 * @param parameter_half_length_z half side of the rectangle in z
 *
 * @return is_outside if out of boundaries
 */
template <size_t TDimension>
bool is_inside_hyper_rectangle(const vec3D_t<int> &point,
                               const size_t &parameter_half_length_x,
                               const size_t &parameter_half_length_y,
                               const size_t &parameter_half_length_z) {
    if constexpr (TDimension == 3) {
        if (abs(point.x) > parameter_half_length_x ||
            abs(point.y) > parameter_half_length_y ||
            abs(point.z) > parameter_half_length_z) {
            return is_outside;
        }
    } else if constexpr (TDimension == 2) {
        if (abs(point.x) > parameter_half_length_x ||
            abs(point.y) > parameter_half_length_y) {
            return is_outside;
        }
    } else if constexpr (TDimension == 1) {
        if (abs(point.x) > parameter_half_length_x) {
            return is_outside;
        }
    } else {
        static_assert(TDimension <= 3, "hyper_rectangle_boundary does not work "
                                       "with dimension higher than 3");
    }
    return is_inside;
}

/**
 * Check if point is inside cube.
 * Assume hyper_cube is centered in 0,0,0
 *
 * @tparam TDimension dimension of the hyper rectangle. options: 1, 2 or 3
 * @param point input point to check if in or out boundaries
 * @param parameter_half_length half side of the cube
 *
 * @return is_outside if out of boundaries
 */
template <size_t TDimension>
bool is_inside_hyper_cube(const vec3D_t<int> &point,
                          const size_t &parameter_half_length) {
    return is_inside_hyper_rectangle<TDimension>(point, parameter_half_length,
                                                 parameter_half_length,
                                                 parameter_half_length);
}

// TODO:Given the normal vector of the 4 planes of the pyramid, use
// convexity to rule out if it is inside or outside.
// template <size_t TDimension>
// bool is_inside_hyper_pyramid(const vec3D_t<int> &point) {
//     return is_inside;
// }

} // namespace boundary
} // namespace perm
#endif
