# Copyright (C) 2019 Pablo Hernandez-Cerdan
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

import _perm as perm
import plot

import unittest

import numpy as np

class TestPerm(unittest.TestCase):
    def test_vec3D(self):
        print("test_vec3D")
        vec0 = perm.vec3Di()
        vec1 = perm.vec3Di(1,1,1)
        plus = vec0 + vec1
        self.assertEqual(plus, vec1)
        # operator []
        vec_op = perm.vec3Di(0,1,2)
        self.assertEqual(vec_op[0], 0)
        self.assertEqual(vec_op[1], 1)
        self.assertEqual(vec_op[2], 2)

    def test_vec3D_numpy(self):
        print("test_vec3D_numpy")
        vec0 = perm.vec3Di()
        arr = np.array(vec0)
        self.assertEqual(arr[0], vec0.x)
        self.assertEqual(arr[1], vec0.y)
        self.assertEqual(arr[2], vec0.z)
        vec0.x = 3
        self.assertEqual(vec0.x, 3)

    def test_single_chain(self):
        print("test_single_chain")
        chain = perm.single_chain()
        self.assertEqual(len(chain.points), 0)
        chain.points.append(perm.vec3Di(1,2,3))
        chain.monomers = 1
        self.assertEqual(len(chain.points), 1)
        self.assertEqual(chain.points[0].x, 1)
        self.assertEqual(chain.points[0].y, 2)
        self.assertEqual(chain.points[0].z, 3)

    def test_random_walk_lattice_2D(self):
        print("test_random_walk_lattice_2D")
        monomers = 5
        chain = perm.random_walk_lattice(N=monomers, dim=2, neighbors=4)
        self.assertEqual(chain.monomers, monomers)
        chain = perm.random_walk_lattice(N=monomers, dim=2, neighbors=8)
        self.assertEqual(chain.monomers, monomers)
        self.assertRaises(RuntimeError, perm.random_walk_lattice, N=monomers, dim=2, neighbors=99)

    def test_random_walk_lattice_3D(self):
        print("test_random_walk_lattice_3D")
        monomers = 5
        chain = perm.random_walk_lattice(N=monomers, dim=3, neighbors=6)
        self.assertEqual(chain.monomers, monomers)
        chain = perm.random_walk_lattice(N=monomers, dim=3, neighbors=18)
        self.assertEqual(chain.monomers, monomers)
        chain = perm.random_walk_lattice(N=monomers, dim=3, neighbors=26)
        self.assertEqual(chain.monomers, monomers)
        self.assertRaises(RuntimeError, perm.random_walk_lattice, N=monomers, dim=3, neighbors=99)

    def test_end_to_end(self):
        print("test_end_to_end")
        chain = perm.single_chain()
        chain.points.append(perm.vec3Di(0,0,0))
        chain.points.append(perm.vec3Di(4,0,3))
        chain.monomers = 2;
        self.assertEqual(chain.ete_vector()[0], 4)
        self.assertEqual(chain.ete_vector()[1], 0)
        self.assertEqual(chain.ete_vector()[2], 3)
        self.assertAlmostEqual(chain.ete_distance(), 5)

    def test_generate_chains_perm(self):
        print("test_generate_chains_perm")
        num_monomers = 20
        parameters_in = perm.parameters_in_t(num_monomers)
        parameters_in.lattice = perm.lattice.d3_n6
        num_chains = 5
        method = "perm"
        parameters_out = perm.generate_chains(num_chains=num_chains, parameters_in=parameters_in, method=method)
        # plot.plot_parameters_out_many_chains(parameters_out)

    def test_generate_chains_rosenbluth(self):
        print("test_generate_chains_rosenbluth")
        num_monomers = 20
        parameters_in = perm.parameters_in_t(num_monomers)
        parameters_in.lattice = perm.lattice.d3_n6
        num_chains = 5
        method = "rosenbluth"
        parameters_out = perm.generate_chains(num_chains=num_chains, parameters_in=parameters_in, method=method)
        # plot.plot_parameters_out_many_chains(parameters_out)

class TestPermPlot(unittest.TestCase):
    def test_plot_chain_2D(self):
        print("plot_chain_2D")
        chain = perm.random_walk_lattice(N=200, dim=2, neighbors=4)
        plot.plot_chain_2D(chain, x=0, y=1)

    def test_plot_chain(self):
        print("plot_chain")
        chain = perm.random_walk_lattice(N=200, dim=2, neighbors=4)
        plot.plot_chain(chain)


if __name__ == "__main__":
    unittest.main()
