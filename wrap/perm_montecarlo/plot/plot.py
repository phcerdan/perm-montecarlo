# Copyright (C) 2019 Pablo Hernandez-Cerdan
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# from .. import perm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_chain_2D(chain, x=0, y=1, ax = "new" ):
    xdata = []
    ydata = []
    for p in chain.points:
        xdata.append(p[x])
        ydata.append(p[y])
    if ax == "new":
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.plot(xdata, ydata)

def plot_chain(chain, ax = "new"):
    xdata = []
    ydata = []
    zdata = []
    for p in chain.points:
        xdata.append(p[0])
        ydata.append(p[1])
        zdata.append(p[2])
    if ax == "new":
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    ax.plot(xdata, ydata, zdata)

def plot_parameters_out_many_chains(parameters_out_many_chains, ax = "new" ):
    if ax == "new":
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    for chain in parameters_out_many_chains.chains:
        plot_chain(chain, ax)
