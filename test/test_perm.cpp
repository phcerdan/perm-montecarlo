/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "perm.hpp"
#include "gmock/gmock.h"

TEST(PERM, random_walk_lattice_2D) {
    size_t monomers = 100;
    const size_t dimension = 2;
    const size_t neighbors = 4;
    auto chain = perm::random_walk_lattice(monomers, dimension, neighbors);
    chain.print(std::cout);
}
