/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_ENERGY_FUNCTIONS_HPP
#define PERM_ENERGY_FUNCTIONS_HPP

#include "perm_common_types.hpp"
#include "single_chain.hpp"
#include "vec3D.hpp"
#include <functional>

namespace perm {

using energy_grow_func_t = std::function<perm::float_t(
        const single_chain_t<int> &, const vec3D_t<int> &)>;

/**
 * Compute total energy of a chain from a energy_grow function.
 *
 * @param input_chain
 * @param energy_func any energy_grow_func_t
 *
 * @return total energy of the input chain
 */
perm::float_t chain_total_energy(const single_chain_t<int> &input_chain,
                                 energy_grow_func_t &energy_func);

perm::float_t energy_grow_zero(const single_chain_t<int> &input_chain,
                               const vec3D_t<int> &new_monomer);
/**
 * Bending energy between last two beads of the chain and new_monomer
 * \f$ \epsilon = k*sum_{j=1}^{N-1}(1-cos_{i,j}) \f$
 *
 * k is related to persistence length lp through:
 *
 * From:
 * Modeling the stretching of wormlike chain in the presence of excluded volume
 * Li et al. Soft Matter 2015
 *
 * \f[
 * \frac{l_p}{w} = \frac{k}{k + 1 - k\coth(k)}
 * \f]
 * where w is the width of the chain
 *
 * @param input_chain
 * @param new_monomer
 * @param k bending energy parameter, related to lp through
 *
 *
 * @return
 */
perm::float_t energy_grow_bending(const single_chain_t<int> &input_chain,
                                  const vec3D_t<int> &new_monomer,
                                  const perm::float_t &k);

} // end namespace perm

#endif
