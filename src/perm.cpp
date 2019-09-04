/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "perm.hpp"
#include "perm_rng.hpp"
#include <iostream>
#include <unordered_set>

namespace perm {

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

void parameters_out_t::print(std::ostream &os) const {
    in.print(os);
    os << "energy= " << energy << std::endl;
    chain.print(os);
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
    while (chain.monomers < parameters_in.monomers) {
        const auto valid_directions = atmosphere_valid_directions(
                chain.points.back(), occupied_map, parameters_in.lattice);
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
            // Choose a direction to grow from the valid atmosphere.
            std::uniform_int_distribution<int> uid_valid(0, num_valid_dirs - 1);
            const auto valid_dir_index = uid_valid(RNG::engine());
            const auto dir = valid_directions[valid_dir_index];
            const auto new_monomer = perm::plus(chain.points.back(),
                                                parameters_in.lattice.at(dir));
            // add monomer
            chain.points.emplace_back(new_monomer);
            chain.monomers++;
            occupied_map.insert(new_monomer);
            weight *= num_valid_dirs;
            // population control
            if (weight > parameters_in.weight_threshold_high[chain.monomers]) {
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

            } else if (weight <
                       parameters_in.weight_threshold_low[chain.monomers]) {
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
        if (chain.monomers < parameters_in.monomers) {
            perm_grow(chain, weight, occupied_map, chain_stack, parameters_in);
        }
        if (chain.monomers == parameters_in.monomers) {
            return {chain, weight};
        }
    }

    // saw not found return empty chain
    return {single_chain_t<int>(), 0};
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

float_t
energy(const single_chain_t<int> &chain,
       const std::function<float_t(const vec3D_t<int> &, const vec3D_t<int> &)>
               &energy_pair_func) {

    if (chain.monomers == 0 || chain.points.empty()) {
        return 0.0;
    }

    auto &points = chain.points;
    double sum = 0.0;
    for (size_t i = 0; i != points.size() - 1; i++) {
        sum += energy_pair_func(points[i], points[i + 1]);
    }
    return sum;
}
} // namespace perm
