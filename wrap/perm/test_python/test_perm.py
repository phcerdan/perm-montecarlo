# Copyright (C) 2019 Pablo Hernandez-Cerdan
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

import _perm as perm
import unittest
import matplotlib.pyplot as plt

import numpy as np

class TestPerm(unittest.TestCase):
    def test_vec3D(self):
        print("test_vec3D")
        vec0 = perm.vec3D()
        print(vec0)
        vec1 = perm.vec3D(1,1,1)
        print(vec1)
        plus = vec0 + vec1
        self.assertEqual(plus, vec1)

    def test_vec3D_numpy(self):
        print("test_vec3D_numpy")
        vec0 = perm.vec3D()
        arr = np.array(vec0)
        print(arr.shape)
        self.assertEqual(arr[0], vec0.x)
        self.assertEqual(arr[1], vec0.y)
        self.assertEqual(arr[2], vec0.z)
        vec0.x = 3
        self.assertEqual(vec0.x, 3)

    def test_single_chain(self):
        print("test_single_chain")
        chain = perm.single_chain()
        self.assertEqual(len(chain.points), 0)
        chain.points.append(perm.vec3D(1,2,3))
        chain.monomers = 1
        self.assertEqual(len(chain.points), 1)
        self.assertEqual(chain.points[0].x, 1)
        self.assertEqual(chain.points[0].y, 2)
        self.assertEqual(chain.points[0].z, 3)
        # operator[]
        self.assertEqual(chain.points[0][0], 1)
        self.assertEqual(chain.points[0][1], 2)
        self.assertEqual(chain.points[0][2], 3)


    def test_random_walk_lattice_2D(self):
        print("test_random_walk_lattice_2D")
        chain = perm.random_walk_lattice_2D(5)
        print(chain)

    def plot_chain(self):
        chain = perm.random_walk_lattice_2D(200)
        xdata = []
        ydata = []
        zdata = []
        for p in chain.points:
            xdata.append(p.x)
            ydata.append(p.y)
            zdata.append(p.z)
        plt.plot(xdata, ydata)


if __name__ == "__main__":
    unittest.main()
