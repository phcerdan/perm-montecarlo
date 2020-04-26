/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "perm.hpp"
#include "perm_rng.hpp"
#include "single_chain.hpp"
#include <algorithm> // for std::find
#include <cmath>     // for exp
#include <iostream>
#include <numeric> // for std::accumulate
#include <string> // for std::to_string
#include <unordered_set>

namespace perm {

bool parameters_in_t::continue_condition(
        const single_chain_t<int> &chain) const {
    return this->monomers != 0 ? chain.monomers < this->monomers
                               : perm::end_to_end_distance(chain) <
                                         this->end_to_end_distance;
}
bool parameters_in_t::done_condition(const single_chain_t<int> &chain) const {
    return this->monomers != 0 ? chain.monomers == this->monomers
                               : perm::end_to_end_distance(chain) >
                                         this->end_to_end_distance;
}
void parameters_in_t::print(std::ostream &os) const {
    os << "max_tries= " << max_tries << std::endl;
    os << "monomers= " << monomers << std::endl;
    os << "weight_threshold_high = " << std::endl;
    for (const auto &wth : weight_threshold_high) {
        os << wth << ", ";
    }
    os << std::endl;
    os << "weight_threshold_low = " << std::endl;
    for (const auto &wtl : weight_threshold_low) {
        os << wtl << ", ";
    }
    os << std::endl;
    os << "end_to_end_distance= " << end_to_end_distance << std::endl;
}

void parameters_out_many_chains_t::print(std::ostream &os) const {
    // in.print(os);
    os << "num_chains: " << num_chains << " (" << chains.size() << ")"
       << std::endl;
    os << "num_weights: " << weights.size() << std::endl;

    const auto N = num_chains;
    const auto average_weight = std::accumulate(
            std::begin(weights), std::end(weights), 0.0,
            [&N](const double &a, const double &b) { return a / N + b / N; });
    os << "average_weight: " << average_weight << std::endl;
}

single_chain_t<int> random_walk_lattice(
        const size_t &monomers,
        const std::function<vec3D_t<int>(void)> &rand_lattice_func) {
    single_chain_t<int> chain;
    chain.points.emplace_back(vec3D_t<int>{0, 0, 0});
    chain.monomers++;
    while (chain.monomers < monomers) {
        chain.points.emplace_back(
                perm::plus(chain.points.back(), rand_lattice_func()));
        chain.monomers++;
    }
    return chain;
}
single_chain_t<int> random_walk_lattice(const size_t &monomers,
                                        const size_t &dimension,
                                        const size_t &neighbors) {
    if (dimension == 1) {
        if (neighbors == 2) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_1D_2n);
        } else {
            throw std::logic_error(
                    "In 1D, neighbors only be 2, but neighbors= " +
                    std::to_string(neighbors));
        }
    } else if (dimension == 2) {
        if (neighbors == 4) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_2D_4n);
        } else if (neighbors == 8) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_2D_8n);
        } else {
            throw std::logic_error(
                    "In 2D, neighbors can be 4 or 8 but neighbors= " +
                    std::to_string(neighbors));
        }
    } else if (dimension == 3) {
        if (neighbors == 6) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_3D_6n);
        } else if (neighbors == 18) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_3D_18n);
        } else if (neighbors == 26) {
            return random_walk_lattice(monomers, &RNG::rand_lattice_3D_26n);
        } else {
            throw std::logic_error(
                    "In 3D, neighbors can be 6, 18 or 26 but neighbors= " +
                    std::to_string(neighbors));
        }
    } else {
        throw std::logic_error(
                "Only dimension 1, 2 or 3 supported but dimension= " +
                std::to_string(dimension));
    }
}

/**
 * Simple sampling discards the whole chain when there is an intersection,
 * starting again from scratch.
 * The chances of getting a SAW decay exponentially with the number of
 * monomers, which makes it inneficient for long chains.
 *
 * @param parameters_in
 *
 * @return
 */
single_chain_t<int> mc_saw_simple_sampling(const size_t &monomers,
                                           const size_t &mc_max_tries,
                                           const lattice_map_t &lattice) {
    auto final_chain = single_chain_t<int>();
    // start at 0
    auto zero_vec = vec3D_t<int>{0, 0, 0};
    final_chain.points.emplace_back(zero_vec);
    final_chain.monomers++;
    auto inverse_direction_lattice =
            perm::get_inverse_direction_lattice(lattice);
    using occupied_map_t = std::unordered_set<vec3D_t<int>, vec3D_int_hasher>;
    static thread_local std::uniform_int_distribution<int> uid(
            0, lattice.size() - 1);
    size_t s = 0;
    for (; s < mc_max_tries; s++) {
        auto occupied_map = occupied_map_t();
        occupied_map.insert(zero_vec);
        // copy the chain
        auto chain = final_chain;
        int last_dir = -100; // init last direction of the bond
        int dir = -100;      // init current direction
        // grow the first edge choosing depending on the lattice
        if (chain.monomers == 1 && monomers > 1) {
            dir = uid(RNG::engine());
            const auto new_monomer =
                    perm::plus(chain.points.back(), lattice.at(dir));
            chain.points.emplace_back(new_monomer);
            chain.monomers++;
            occupied_map.insert(new_monomer);
            last_dir = dir;
        }
        bool discard_chain = false;
        while (chain.monomers < monomers) {
            // Choose a direction to grow, different than the inverse of the
            // last movement (to speed-up, and not discard the whole chain just
            // yet)
            do {
                dir = uid(RNG::engine());
            } while (dir == inverse_direction_lattice.at(last_dir));
            last_dir = dir;
            auto new_monomer = perm::plus(chain.points.back(), lattice.at(dir));
            // Check that the new monomer does not collide with any other
            auto search = occupied_map.find(new_monomer);
            if (search == occupied_map.end()) {
                chain.points.emplace_back(new_monomer);
                chain.monomers++;
                occupied_map.insert(new_monomer);
            } else {
                // discard chain and starts again
                discard_chain = true;
                break;
            }
        }
        if (!discard_chain) {
            // You got a saw chain with the right amount of monomers
            // stop trying and return the chain.
            final_chain = chain;
            break;
        }
    }

    if (s == mc_max_tries) {
        // saw not found, return empty chain
        return single_chain_t<int>();
    } else {
        return final_chain;
    }
} // namespace perm

std::vector<int> atmosphere_valid_directions(
        const vec3D_t<int> &monomer,
        const std::unordered_set<vec3D_t<int>, vec3D_int_hasher> &occupied_map,
        const lattice_map_t &lattice) {
    std::vector<int> valid_directions;
    for (const auto &key_vector_pair : lattice) {
        const auto &key = key_vector_pair.first;
        const auto &dir_vector = key_vector_pair.second;
        const auto new_monomer = perm::plus(monomer, dir_vector);
        const auto search = occupied_map.find(new_monomer);
        if (search == occupied_map.end()) { // not occupied
            valid_directions.emplace_back(key);
        }
    }
    return valid_directions;
}

std::vector<int> valid_directions_with_atmosphere_and_bondary_condition(
        const vec3D_t<int> &monomer,
        const std::unordered_set<vec3D_t<int>, vec3D_int_hasher> &occupied_map,
        const lattice_map_t &lattice,
        const boundary_func_t &is_inside_boundary_func) {
    auto valid_directions =
            atmosphere_valid_directions(monomer, occupied_map, lattice);
    valid_directions.erase(
            std::remove_if(valid_directions.begin(), valid_directions.end(),
                           [&is_inside_boundary_func, &monomer,
                            &lattice](const int &dir) {
                               return !is_inside_boundary_func(
                                       perm::plus(monomer, lattice.at(dir)));
                           }),
            valid_directions.end());
    return valid_directions;
}

size_t non_bonded_nearest_neighbors(
        const std::vector<vec3D_t<int>> &chain_points,
        const vec3D_t<int> &monomer,
        const std::unordered_set<vec3D_t<int>, vec3D_int_hasher> &occupied_map,
        const lattice_map_t &lattice) {

    const std::vector<int> valid_directions =
            atmosphere_valid_directions(monomer, occupied_map, lattice);
    const size_t neighbor_monomers = valid_directions.size();
    const size_t lattice_size = lattice.size();
    const size_t touching_pairs = lattice_size - neighbor_monomers;
    const auto search = std::find(std::begin(chain_points),
                                  std::end(chain_points), monomer);
    if (search == std::end(chain_points)) { // monomer not in the chain
        return touching_pairs;              // all touching pairs are non-bonded
    } else { // monomer is in the chain, (guarantees non-empty chain)
        if (monomer == chain_points.back() || monomer == chain_points[0]) {
            // only one bond, one pair is bonded
            return touching_pairs - 1;
        } else {
            // monomer is in the middle of the chain, two pairs are bonded
            return touching_pairs - 2;
        }
    }
}

/**
 * Rosenbluth sampling
 * The new monomer is selected from a set of nearest neighbors that are not
 * occupied already (these neighbors are in the endpoint atmosphere of the walk)
 *
 * It returns a chain and a weight associated to the chain.
 * where the weight are the products of valid movements for each monomer.
 * \f[
 * W = \prod_{i=0}^{n-1}W_i = \prod_{i=0}^{n-1}\sigma_i
 * \f]
 *
 * where \f$n\f$ is the number of monomers \f$\sigma_i = a_{+}^e\f$ is the
 * positive atmosphere of the monomer i, i.e, the number of open/valid
 * future movements.
 *
 * @param monomers
 * @param mc_max_tries
 *
 * @return single_chain_t, and weight
 */
std::pair<single_chain_t<int>, double>
mc_saw_rosenbluth_sampling(const size_t &monomers,
                           const size_t &mc_max_tries,
                           const lattice_map_t &lattice) {
    auto final_chain = single_chain_t<int>();
    double final_weight = 1;
    // start at 0
    auto zero_vec = vec3D_t<int>{0, 0, 0};
    final_chain.points.emplace_back(zero_vec);
    final_chain.monomers++;
    using occupied_map_t = std::unordered_set<vec3D_t<int>, vec3D_int_hasher>;
    static thread_local std::uniform_int_distribution<int> uid(
            0, lattice.size() - 1);
    size_t s = 0;
    for (; s < mc_max_tries; s++) {
        auto occupied_map = occupied_map_t();
        occupied_map.insert(zero_vec);
        // copy the chain
        auto chain = final_chain;
        auto weight = final_weight;
        int dir = -100; // current direction
        // grow the first edge choosing depending on the lattice
        if (chain.monomers == 1 && monomers > 1) {
            dir = uid(RNG::engine());
            const auto new_monomer =
                    perm::plus(chain.points.back(), lattice.at(dir));
            chain.points.emplace_back(new_monomer);
            chain.monomers++;
            occupied_map.insert(new_monomer);
            weight *= lattice.size();
        }
        bool discard_chain = false;
        while (chain.monomers < monomers) {
            // Compute the number of valid movements (i.e. atmosphere of
            // endpoint)
            const auto valid_directions = atmosphere_valid_directions(
                    chain.points.back(), occupied_map, lattice);
            const auto num_valid_dirs = valid_directions.size();
            if (valid_directions.empty()) {
                // discard chain and starts again
                discard_chain = true;
                break;
            }
            // Choose a direction to grow from the valid atmosphere.
            std::uniform_int_distribution<int> uid_valid(0, num_valid_dirs - 1);
            const auto valid_dir_index = uid_valid(RNG::engine());
            const auto dir = valid_directions[valid_dir_index];
            const auto new_monomer =
                    perm::plus(chain.points.back(), lattice.at(dir));
            chain.points.emplace_back(new_monomer);
            chain.monomers++;
            occupied_map.insert(new_monomer);
            weight *= num_valid_dirs;
        }

        if (!discard_chain) {
            // You got a saw chain with the right amount of monomers
            // stop trying and return the chain.
            final_chain = chain;
            final_weight = weight;
            break;
        }
    }

    if (s == mc_max_tries) {
        // saw not found, return empty chain
        return {single_chain_t<int>(), 0};
    } else {
        return {final_chain, final_weight};
    }
}

// Initialize references outside
void perm_grow(single_chain_t<int> &chain,
               double &weight,
               occupied_map_t &occupied_map,
               chain_stack_t &chain_stack,
               const parameters_in_t &parameters_in) {
    while (parameters_in.continue_condition(chain)) {
        // Check if monomer can be added at the end, or all possible sites
        // are occupied by other monomers
        const auto valid_directions =
                valid_directions_with_atmosphere_and_bondary_condition(
                        chain.points.back(), occupied_map,
                        parameters_in.lattice,
                        parameters_in.is_inside_boundary_func);
        // const auto valid_directions = atmosphere_valid_directions(
        //         chain.points.back(), occupied_map, parameters_in.lattice);
        // Check if any go out of the boundaries of the volume
        const auto num_valid_dirs = valid_directions.size();
        if (valid_directions.empty()) {
            // discard chain, keep growing next chain from the chain_stack
            // or abort this tour if stack is empty
            if (chain_stack.size() > 1) {
                chain_stack.pop();
                auto &top = chain_stack.top();
                auto &chain = top.chain;
                auto &weight = top.weight;
                auto &occupied_map = top.occupied_map;
                perm_grow(chain, weight, occupied_map, chain_stack,
                          parameters_in);
            } else {
                break;
            }
        } else { // this chain can still grow
            // // Choose a direction to grow from the valid atmosphere.
            // std::uniform_int_distribution<int> uid_valid(0, num_valid_dirs -
            // 1); const auto valid_dir_index = uid_valid(RNG::engine()); const
            // auto dir = valid_directions[valid_dir_index]; const auto
            // new_monomer = perm::plus(chain.points.back(),
            //                                     parameters_in.lattice.at(dir));

            // Choose a valid direction, but weighted with the bending energy
            // To recover a non-weighted, use energy_grow_zero as function.
            std::vector<perm::float_t> energies(num_valid_dirs);
            std::vector<perm::float_t> weights_choose_dir(num_valid_dirs);
            for (size_t i = 0; i < num_valid_dirs; i++) {
                const auto new_valid_monomer = perm::plus(
                        chain.points.back(),
                        parameters_in.lattice.at(valid_directions[i]));
                energies[i] = parameters_in.energy_grow_func(chain,
                                                             new_valid_monomer);
                weights_choose_dir[i] =
                        1.0 / exp(parameters_in.beta * energies[i]);
            }
            std::discrete_distribution<int> uid_weighted(
                    std::begin(weights_choose_dir),
                    std::end(weights_choose_dir));
            const auto valid_dir_index = uid_weighted(RNG::engine());
            const auto dir = valid_directions[valid_dir_index];
            const auto new_monomer = perm::plus(chain.points.back(),
                                                parameters_in.lattice.at(dir));

            // From: A review of Monte Carlo simulations of polymers with PERM
            // 2011, Hsu, Grassberg
            // non-bonded interaction: q = exp(-\beta*\epsilon)
            // partition sum: Z_N(q) = \sum_walks q^m
            // m = total number of non-bonded interactions
            // w_n = q^{m_n} / p_n (9)
            // p_n = 1 / n_{free}
            // n_free = free positions for the N-1 monomer
            // m_n = number of neighbors of the new site already occupied for
            // non-bonded monomers. m = \sum_{n=0}^{N} m_n
            // W_N = w_N W_{N-1} = \prod_{n=0}^{N} w_n
            // perm::float_t q = exp(beta * parameters_in.energy_grow_func(
            //         chain, new_monomer)); // eq (9)
            // perm::float_t p_n_bias = 1.0/num_valid_dirs // eq (9)
            // perm::float_t weight_n = std::pow(q, m_n) / p_n_bias // eq (9)
            weight *= num_valid_dirs *
                      exp(parameters_in.beta *
                          parameters_in.energy_grow_func(chain, new_monomer));
            chain.points.emplace_back(new_monomer);
            chain.monomers++;
            occupied_map.insert(new_monomer);
            // population control
            // Check for existence of the input weights first
            const bool weight_threshold_high_monomer_exists =
                    chain.monomers < parameters_in.weight_threshold_high.size();
            const bool weight_threshold_low_monomer_exists =
                    chain.monomers < parameters_in.weight_threshold_low.size();
            if (weight_threshold_high_monomer_exists &&
                weight > parameters_in.weight_threshold_high[chain.monomers]) {
                // enrichment
                // TODO num_copies might be a external parameter
                const size_t num_copies = 1;
                const double weight_reduction = 1.0 / (num_copies + 1);
                chain_info_t chain_info;
                chain_info.chain = chain;
                chain_info.occupied_map = occupied_map;
                // the weight of the new chain is reduced
                chain_info.weight = weight * weight_reduction;
                chain_stack.push(chain_info);
                // and also reduce the weight of the current chain
                weight *= weight_reduction;

            } else if (weight_threshold_low_monomer_exists &&
                       weight < parameters_in
                                        .weight_threshold_low[chain.monomers]) {
                // prune
                // TODO probability of kill might be a external parameter
                const double kill_probability = 0.5;
                const double weight_increase = 2; // 1.0 / kill_probability
                const bool kill_configuration =
                        RNG::random_bool(kill_probability);
                if (kill_configuration) {
                    if (chain_stack.size() > 1) {
                        chain_stack.pop();
                        auto &top = chain_stack.top();
                        auto &chain = top.chain;
                        auto &weight = top.weight;
                        auto &occupied_map = top.occupied_map;
                        perm_grow(chain, weight, occupied_map, chain_stack,
                                  parameters_in);
                    } else {
                        break;
                    }
                } else {
                    // keep chain but increase its weight
                    weight *= weight_increase;
                    perm_grow(chain, weight, occupied_map, chain_stack,
                              parameters_in);
                }
            } else { // keep growing without variation
                // std::cout << chain.monomers << std::endl;
                perm_grow(chain, weight, occupied_map, chain_stack,
                          parameters_in);
            }
        }
    }
}

std::pair<single_chain_t<int>, double>
mc_saw_perm(const parameters_in_t &parameters_in) {
    using occupied_map_t = std::unordered_set<vec3D_t<int>, vec3D_int_hasher>;
    size_t s = 0;
    for (; s < parameters_in.max_tries; s++) {
        // std::cout << "try: s:  " << s << std::endl;
        chain_stack_t chain_stack;
        chain_info_t chain_info;
        chain_info.chain = single_chain_t<int>();
        chain_info.weight = 1;
        chain_info.occupied_map = occupied_map_t();
        chain_stack.push(chain_info);
        auto &top = chain_stack.top();
        auto &chain = top.chain;
        double &weight = top.weight;
        auto &occupied_map = top.occupied_map;
        // start at 0
        auto zero_vec = vec3D_t<int>{0, 0, 0};
        chain.points.emplace_back(zero_vec);
        chain.monomers++;
        occupied_map.insert(zero_vec);
        if (parameters_in.continue_condition(chain)) {
            perm_grow(chain, weight, occupied_map, chain_stack, parameters_in);
        }
        if (parameters_in.done_condition(chain)) {
            return {chain, weight};
        }
    }

    // saw not found return empty chain
    return {single_chain_t<int>(), 0};
}

} // namespace perm
