/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_COMMON_PY_HPP
#define PERM_COMMON_PY_HPP

#include "perm_common_types.hpp"
#include <pybind11/stl_bind.h>
PYBIND11_MAKE_OPAQUE(std::vector<perm::vec3D_t<int>>);
PYBIND11_MAKE_OPAQUE(std::vector<perm::float_t>);
#endif
