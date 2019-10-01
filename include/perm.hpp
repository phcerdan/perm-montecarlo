/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PERM_HPP
#define PERM_HPP

#include "energy_functions.hpp"
#include "lattice_boundary.hpp" // for boundary_func_t
#include "lattice_lut.hpp"      // for lattice_map_t
#include "perm_common_types.hpp"
#include "single_chain.hpp"
#include "vec3D.hpp"
#include <cmath>
#include <cstddef> // For size_t
#include <functional>
#include <limits>
#include <ostream>
#include <stack>         // for chain_stack_t
#include <unordered_set> // for occupied_map_t

namespace perm {
struct parameters_in_t {

    size_t max_tries = 1000;
    size_t monomers = 100;
    std::vector<perm::float_t> weight_threshold_high =
            std::vector<perm::float_t>(
                    monomers, std::numeric_limits<perm::float_t>::infinity());
    std::vector<perm::float_t> weight_threshold_low =
            std::vector<perm::float_t>(monomers, 0.0);
    lattice_map_t lattice = lattice_2D_4n;
    /// Used only if monomers = 0
    float_t end_to_end_distance = 0.0;
    bool continue_condition(const single_chain_t<int> &) const;
    bool done_condition(const single_chain_t<int> &) const;
    void print(std::ostream &os) const;
    /**
     * energy function for adding a monomer to an existing chain
     *
     * @param chain
     * @param new_monomer
     *
     * @return
     */
    energy_grow_func_t energy_grow_func =
            [](const single_chain_t<int> &,
               const vec3D_t<int> &) -> perm::float_t { return 0.0; };
    /// \f$ kT = \frac{1}{\beta} \f$
    float_t kT = 1.0;
    float_t beta = 1.0;

    /**
     * Bondary_function to define the enclosing volume
     * Can be externally defined by user with a lambda function, for example:
     * parameters_in.boundary_func =
     *  [&parameter_half_length = parameter_half_length](
     *       const vec3D_t<int> & point) -> bool {
     *          return is_inside_hyper_cube(point, parameter_half_length);
     * };
     *
     */
    boundary_func_t is_inside_boundary_func = [](const vec3D_t<int> &) -> bool {
        return perm::boundary::is_inside;
    };

    // Constructors/Destructors
    explicit parameters_in_t(const size_t &num_monomers) : parameters_in_t() {
        monomers = num_monomers;
        weight_threshold_high = std::vector<perm::float_t>(
                num_monomers, std::numeric_limits<perm::float_t>::infinity());
        weight_threshold_low = std::vector<perm::float_t>(num_monomers, 0.0);
    };
    parameters_in_t() = default;
    parameters_in_t(const parameters_in_t &) = default;
    parameters_in_t(parameters_in_t &&) = default;
    parameters_in_t &operator=(const parameters_in_t &) = default;
    parameters_in_t &operator=(parameters_in_t &&) = default;
    virtual ~parameters_in_t() = default;
};

struct parameters_out_many_chains_t {
    // parameters_out_many_chains_t() = default;
    // parameters_out_many_chains_t(const parameters_in_t &input) : in(input){};
    // parameters_in_t in;
    size_t num_chains;
    std::vector<single_chain_t<int>> chains;
    std::vector<double> weights;
    void print(std::ostream &os) const;
};

single_chain_t<int>
random_walk_lattice(const size_t &monomers,
                    const std::function<vec3D_t<int>(void)> &rand_lattice_func);
single_chain_t<int> random_walk_lattice(const size_t &monomers,
                                        const size_t &dimension,
                                        const size_t &neighbors);

single_chain_t<int> mc_saw_simple_sampling(const size_t &monomers,
                                           const size_t &mc_max_tries,
                                           const lattice_map_t &lattice);
/**
 * Compute the directions of growth that are valid for the input monomer.
 *
 * @param monomer input vector
 * @param occupied_map occupied_map is the set of already occupied places.
 * @param lattice lattice maps are already defined, they map
 * integers with vector directions, for example lattice_2D_4n;
 *
 * @return vector with valid directions (integers), use lattice again to recover
 * the vector.
 */
std::vector<int> atmosphere_valid_directions(
        const vec3D_t<int> &monomer,
        const std::unordered_set<vec3D_t<int>, vec3D_int_hasher> &occupied_map,
        const lattice_map_t &lattice);

std::vector<int> valid_directions_with_atmosphere_and_bondary_condition(
        const vec3D_t<int> &monomer,
        const std::unordered_set<vec3D_t<int>, vec3D_int_hasher> &occupied_map,
        const lattice_map_t &lattice,
        const boundary_func_t &is_inside_boundary_func);

std::pair<single_chain_t<int>, double>
mc_saw_rosenbluth_sampling(const size_t &monomers,
                           const size_t &mc_max_tries,
                           const lattice_map_t &lattice);

using occupied_map_t = std::unordered_set<vec3D_t<int>, vec3D_int_hasher>;
struct chain_info_t {
    single_chain_t<int> chain;
    double weight;
    occupied_map_t occupied_map;
};
// stack only allows to push/pop from the top, and references are always valid
using chain_stack_t = std::stack<chain_info_t>;

/**
 * Recursive function with perm algorithm. Use @mc_saw_perm
 *
 * @param chain
 * @param weight
 * @param occupied_map
 * @param chain_stack
 * @param parameters_in
 */
void perm_grow(single_chain_t<int> &chain,
               double &weight,
               occupied_map_t &occupied_map,
               chain_stack_t &chain_stack,
               const parameters_in_t &parameters_in);

std::pair<single_chain_t<int>, double>
mc_saw_perm(const parameters_in_t &parameters_in);

template <typename Func, typename... Args>
parameters_out_many_chains_t
generate_chains(const size_t &num_chains, const Func &func, Args &&... args) {
    parameters_out_many_chains_t out;
    out.num_chains = num_chains;
    out.chains.reserve(num_chains);
    out.weights.reserve(num_chains);
    for (size_t i = 0; i < num_chains; i++) {
        auto [chain, weight] = func(std::forward<Args>(args)...);
        out.chains.emplace_back(chain);
        out.weights.emplace_back(weight);
    }
    return out;
}

} // namespace perm
#endif
