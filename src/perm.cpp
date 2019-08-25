/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#include "perm.hpp"
#include <iostream>

namespace perm {
    using precision_t = float;

    struct parameters {
        size_t steps = 1000;
        size_t num_monomers = 100;
        /// Used only if num_monomers = 0
        precision_t end_to_end_distance = 0.0;
    }

    void hola() {
        std::cout << "HOLA" << std::endl;
    }
}
