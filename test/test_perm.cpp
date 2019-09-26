/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "lattice_boundary.hpp"
#include "lattice_lut.hpp"
#include "perm.hpp"
#include "gmock/gmock.h"
#include <numeric> // for accumulate

TEST(PERM, vec3D) {
    perm::vec3D_t<int> v = {1, 2, 3};
    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 3);
    EXPECT_EQ(v.x, 1);
    EXPECT_EQ(v.y, 2);
    EXPECT_EQ(v.z, 3);
    EXPECT_EQ(v.x - v.x, 0);
}

TEST(PERM, random_walk_lattice_2D) {
    size_t monomers = 5;
    const size_t dimension = 2;
    const size_t neighbors = 4;
    auto chain = perm::random_walk_lattice(monomers, dimension, neighbors);
    chain.print(std::cout);
}

TEST(PERM, random_walk_lattice_has_expected_end_to_end_distance) {
    const size_t monomers = 100;
    const size_t bonds = monomers - 1;
    const size_t dimension = 2;
    const size_t neighbors = 4;
    const size_t num_experiments = 100;
    auto ete_distances = std::vector<double>(num_experiments);
    auto gyration_radiuses_square = std::vector<double>(num_experiments);
    for (int i = 0; i < num_experiments; i++) {
        const auto chain =
                perm::random_walk_lattice(monomers, dimension, neighbors);
        ete_distances[i] = perm::end_to_end_distance(chain);
        gyration_radiuses_square[i] = perm::gyration_radius_square(chain);
    }
    const auto average_ete_distance =
            std::accumulate(ete_distances.begin(), ete_distances.end(), 0.0) /
            ete_distances.size();
    const auto bond_length_lattice =
            perm::bond_length_lattice(dimension, neighbors);
    const auto expected_ete_distance = sqrt(monomers) * bond_length_lattice;
    // XXX why is not closer to the expected theoretical value?
    EXPECT_GT(average_ete_distance / expected_ete_distance, 0.8)
            << " average_ete_distance: " << average_ete_distance
            << "\n expected_ete_distance: " << expected_ete_distance;
    // gyration_2 = N*b^2/6.0
    const auto average_gyration_radius_square =
            std::accumulate(gyration_radiuses_square.begin(),
                            gyration_radiuses_square.end(), 0.0) /
            gyration_radiuses_square.size();
    const auto expected_gyration_radius_square =
            monomers * bond_length_lattice * bond_length_lattice / 6.0;
    EXPECT_GT(average_gyration_radius_square / expected_gyration_radius_square,
              0.8)
            << " average_gyration_radius_square: "
            << average_gyration_radius_square
            << "\n expected_gyration_radius_square: "
            << expected_gyration_radius_square;
}

TEST(PERM, mc_saw_simple_sampling) {
    const size_t monomers = 5;
    const size_t tries = 1000;
    const auto lattice = perm::lattice_2D_4n;
    auto chain = perm::mc_saw_simple_sampling(monomers, tries, lattice);
    EXPECT_FALSE(chain.points.empty());
    chain.print(std::cout);
}

TEST(PERM, mc_saw_rosenbluth_sampling) {
    const size_t monomers = 10;
    const size_t tries = 1000;
    const auto lattice = perm::lattice_2D_4n;
    auto chain_weight_pair =
            perm::mc_saw_rosenbluth_sampling(monomers, tries, lattice);
    auto &chain = chain_weight_pair.first;
    auto &weight = chain_weight_pair.second;
    EXPECT_FALSE(chain.points.empty());
    chain.print(std::cout);
}

TEST(PERM, mc_saw_perm) {
    const size_t monomers = 10;
    const size_t tries = 1000;
    const auto lattice = perm::lattice_2D_4n;
    perm::parameters_in_t parameters_in;
    parameters_in.monomers = monomers;
    parameters_in.max_tries = tries;
    parameters_in.lattice = lattice;
    std::cout << "Parameters = " << std::endl;
    parameters_in.print(std::cout);
    auto chain_weight_pair = perm::mc_saw_perm(parameters_in);
    auto &chain = chain_weight_pair.first;
    auto &weight = chain_weight_pair.second;
    EXPECT_FALSE(chain.points.empty());
    chain.print(std::cout);
    std::cout << "Weight = " << weight << std::endl;

    const size_t num_experiments = 10;
    std::vector<perm::single_chain_t<int>> chains;
    std::vector<perm::float_t> weights;
    for (int i = 0; i < num_experiments; i++) {
        auto chain_weight_pair = perm::mc_saw_perm(parameters_in);
        chains.push_back(chain_weight_pair.first);
        weights.push_back(chain_weight_pair.second);
    }
    EXPECT_EQ(chains.size(), num_experiments);
    EXPECT_EQ(weights.size(), num_experiments);
}

TEST(lattice_boundary, hyper_rectangle) {
    constexpr size_t dimension = 3;
    const size_t rectangle_half_length_x = 5;
    const size_t rectangle_half_length_y = 3;
    const size_t rectangle_half_length_z = 2;
    perm::vec3D_t<int> inside_point{4, 2, 1};
    EXPECT_TRUE(perm::boundary::is_inside_hyper_rectangle<dimension>(
            inside_point, rectangle_half_length_x, rectangle_half_length_y,
            rectangle_half_length_z));
    perm::vec3D_t<int> outside_point{6, 3, 1};
    EXPECT_FALSE(perm::boundary::is_inside_hyper_rectangle<dimension>(
            outside_point, rectangle_half_length_x, rectangle_half_length_y,
            rectangle_half_length_z));
    perm::vec3D_t<int> border_point{5, 2, 1};
    EXPECT_TRUE(perm::boundary::is_inside_hyper_rectangle<dimension>(
            border_point, rectangle_half_length_x, rectangle_half_length_y,
            rectangle_half_length_z));
    // The following does not compile for dimension 4, as expected.
    // perm::boundary::is_inside_hyper_cube<4>(inside_point,
    // rectangle_half_length_x);
}

TEST(generate_chains, perm) {
    const size_t num_chains = 10;
    const size_t num_monomers = 100;
    auto parameters_in = perm::parameters_in_t(num_monomers);
    auto out = generate_chains(num_chains, perm::mc_saw_perm, parameters_in);
    EXPECT_EQ(out.num_chains, num_chains);
    EXPECT_EQ(out.chains[0].points.size(), num_monomers);
    EXPECT_EQ(out.chains.back().points.size(), num_monomers);
    out.print(std::cout);
}

TEST(generate_chains, rosenbluth) {
    const size_t num_chains = 10;
    const size_t num_monomers = 100;
    const size_t max_tries = 10000;
    const auto lattice = perm::lattice_3D_6n;
    auto out = generate_chains(num_chains, perm::mc_saw_rosenbluth_sampling,
            num_monomers, max_tries, lattice);
    EXPECT_EQ(out.num_chains, num_chains);
    EXPECT_EQ(out.chains[0].points.size(), num_monomers);
    EXPECT_EQ(out.chains.back().points.size(), num_monomers);
    out.print(std::cout);
}
