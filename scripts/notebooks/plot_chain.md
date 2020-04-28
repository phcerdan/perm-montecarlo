```python
%matplotlib inline
```


```python
import random
import statistics
import numpy as np
import scipy.stats
import scipy
import perm
# import perm._perm as perm
import perm.plot as plot
import matplotlib.pyplot as plt
import math
```


```python
# generate_chains (test)
num_monomers = 9
parameters_in = perm.parameters_in_t(num_monomers)
parameters_in.end_to_end_distance = 20.0
parameters_in.lattice = perm.lattice.d3_n26
k_bending=8
energy_grow_func = perm.energy_grow_bending(k=k_bending)
parameters_in.energy_grow_func = energy_grow_func
num_chains = 12
method = "perm"
parameters_out = perm.generate_chains(num_chains=num_chains, parameters_in=parameters_in, method=method)
plot.plot_parameters_out_many_chains(parameters_out)

method = "rosenbluth"
parameters_out = perm.generate_chains(num_chains=num_chains, parameters_in=parameters_in, method=method)
plot.plot_parameters_out_many_chains(parameters_out)
print(parameters_out)
```

    num_chains: 12 (12)
    num_weights: 12
    average_weight: 1.26925e+10
    



![png](scripts/notebooks/plot_chain_files/plot_chain_2_1.png)



![png](scripts/notebooks/plot_chain_files/plot_chain_2_2.png)



```python

```


```python
# Draw one chain
# chain = perm.random_walk_lattice(N=20000, dim=3, neighbors=6)
# plot.plot_chain_2D(chain, 0, 1)
# plot.plot_chain_2D(chain, 0, 2)
# plot.plot_chain_2D(chain, 1, 2)
# plot.plot_chain(chain)
# chain = perm.random_walk_lattice(N=20000, dim=2, neighbors=8)
# plot.plot_chain_2D(chain)
```

# Draw many chains


```python
num_experiments = 1000

monomers = 20
bonds = monomers - 1
dim = 2
neighbors = 8
fig = plt.figure()
if dim == 2:
    ax = fig.add_subplot(111)
else:
    ax = fig.add_subplot(111, projection='3d')

ete_distances = []
for i in range(num_experiments):
    chain = perm.random_walk_lattice(N=monomers, dim=dim, neighbors=neighbors)
    ete_distances.append(chain.ete_distance())
    if dim == 2:
        plot.plot_chain_2D(chain, 0, 1, ax)
    else:
        plot.plot_chain(chain, ax)
average_ete_distance = sum(ete_distances)/num_experiments
expected_ete_distance = math.sqrt(bonds) * perm.bond_length_lattice(dim=dim, neighbors=neighbors)
print("average_ete_distance", average_ete_distance)
print("expected_ete_distance", expected_ete_distance)
print("ratio", average_ete_distance/expected_ete_distance)

```

    average_ete_distance 4.770488355006554
    expected_ete_distance 5.261656473254825
    ratio 0.9066514279780722



![png](scripts/notebooks/plot_chain_files/plot_chain_6_1.png)


# STATS


```python

num_experiments = 1000
monomers = 2000
bonds = monomers - 1
dim = 3
neighbors = 6

ete_distances = []
ete_vectors = []
for _ in range(num_experiments):
    chain = perm.random_walk_lattice(N=monomers, dim=dim, neighbors=neighbors)
    ete_distances.append(chain.ete_distance())
    ete_vectors.append(chain.ete_vector())

average_ete_distance = statistics.mean(ete_distances)
expected_ete_distance = math.sqrt(bonds) * perm.bond_length_lattice(dim=dim, neighbors=neighbors)
print("average_ete_distance", average_ete_distance)
print("expected_ete_distance", expected_ete_distance)
print("ratio", average_ete_distance/expected_ete_distance)
x_positions = [d[0] for d in ete_vectors]
y_positions = [d[1] for d in ete_vectors]
z_positions = [d[2] for d in ete_vectors]
average_pos_x = statistics.mean(x_positions)
average_pos_y = statistics.mean(y_positions)
average_pos_z = statistics.mean(z_positions)
print("average_pos", [average_pos_x, average_pos_y, average_pos_z])
print("expected_pos", [0, 0, 0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
nbins = 100
hist_x, bins_x = np.histogram(x_positions, bins=nbins)
hist_y, bins_y = np.histogram(y_positions, bins=nbins)
hist_z, bins_z = np.histogram(z_positions, bins=nbins)
bins_x_centers = (bins_x[:-1] + bins_x[1:])/2
bins_y_centers = (bins_y[:-1] + bins_y[1:])/2
bins_z_centers = (bins_z[:-1] + bins_z[1:])/2

ax.bar(bins_x_centers, hist_x, zs=-1, zdir='y', label='x_pos')
ax.bar(bins_y_centers, hist_y, zs=0, zdir='y', label='y_pos')
ax.bar(bins_z_centers, hist_z, zs=1, zdir='y', label='z_pos')
ax.set_xlabel('displacement')
ax.set_ylabel('---')
ax.set_zlabel('counts')
ax.legend(loc='best')

fig = plt.figure()
ay = fig.add_subplot(111, projection='3d')
distro_x = scipy.stats.rv_histogram([hist_x, bins_x])
distro_y = scipy.stats.rv_histogram([hist_y, bins_y])
distro_z = scipy.stats.rv_histogram([hist_z, bins_z])
ay.plot(bins_x_centers, distro_x.pdf(bins_x_centers), zs=-1, zdir='y', label='xpos')
ay.plot(bins_y_centers, distro_y.pdf(bins_y_centers), zs=0, zdir='y', label='ypos')
ay.plot(bins_z_centers, distro_z.pdf(bins_z_centers), zs=1, zdir='y', label='zpos')
ay.set_xlabel('displacement')
ay.set_ylabel('---')
ay.set_zlabel('pdf (probability)')
ay.legend(loc='best')
```

    average_ete_distance 41.8280957430458
    expected_ete_distance 44.710177812216315
    ratio 0.9355385683931896
    average_pos [-0.195, 0.864, 0.047]
    expected_pos [0, 0, 0]





    <matplotlib.legend.Legend at 0x7f972db75c40>




![png](scripts/notebooks/plot_chain_files/plot_chain_8_2.png)



![png](scripts/notebooks/plot_chain_files/plot_chain_8_3.png)



```python
# # example pure python
```


```python

def bmean(neighbors):
    """Correct average kuhn monomer length for second nearest neighbors lattices"""
    bmean = 1
    if neighbors == 8:  # 2D
        bmean = (4*1 + 4*math.sqrt(2))/neighbors
    elif neighbors == 18:
        bmean = (6*1 + 12*math.sqrt(2))/neighbors
    elif neighbors == 26:
        bmean = (6*1 + 12*math.sqrt(2) + 8*math.sqrt(3))/neighbors

    return bmean

def choice_2d_4n():
    return random.choice(
        [np.array([-1, 0, 0]), np.array([0, -1, 0]),
         np.array([1, 0, 0]), np.array([0, 1, 0])])

def choice_3d_6n():
    return random.choice(
        [np.array([-1, 0, 0]), np.array([0, -1, 0]),
         np.array([1, 0, 0]), np.array([0, 1, 0]),
         np.array([0, 0, -1]), np.array([0, 0, 1])]
    )

def random_walk_lattice(monomers, dim):
    chain = []
    chain.append(np.array([0, 0, 0]))
    current_monomers = 1
    if(dim == 2):
        random_func = choice_2d_4n
    elif (dim == 3):
        random_func = choice_3d_6n
    else:
        raise "wrong dim"

    while (current_monomers < monomers):
        new_pos = random_func()
        chain.append(chain[-1] + new_pos)
        current_monomers += 1
    return chain


num_experiments = 1000
monomers = 200
bonds = monomers - 1
dim = 3
neighbors = 6
ete_distances = []
for _ in range(num_experiments):
    chain = random_walk_lattice(monomers, dim=dim)
    ete_distances.append(np.linalg.norm(chain[-1] - chain[0]))
    #print(chain[-1], ete_distances[-1])
# print(chain)

average_ete_distance = statistics.mean(ete_distances)
expected_ete_distance = math.sqrt(bonds)*bmean(neighbors)
print("average_ete_distance", average_ete_distance)
print("expected_ete_distance", expected_ete_distance)
print("ratio", average_ete_distance/expected_ete_distance)

```

    average_ete_distance 13.109292919296976
    expected_ete_distance 14.106735979665885
    ratio 0.9292931361438486


# SAW: Self Avoiding Walk
### simple sampling


```python
monomers = 15
mc_tries = 1000
lattice = perm.lattice.d2_n4
saw_chain = perm.mc_saw_simple_sampling(monomers=monomers, tries=mc_tries, lattice=lattice)
if saw_chain.monomers == 0:
    print("MC did not found a self avoiding walk.")
else:
    print(saw_chain)

plot.plot_chain_2D(saw_chain)
```

    chain.monomers= 15
    0 0 0
    0 -1 0
    1 -1 0
    1 -2 0
    1 -3 0
    2 -3 0
    3 -3 0
    3 -2 0
    3 -1 0
    3 0 0
    3 1 0
    2 1 0
    1 1 0
    1 0 0
    2 0 0
    



![png](scripts/notebooks/plot_chain_files/plot_chain_12_1.png)


### Rosenbluth sampling


```python
monomers = 16
mc_tries = 1000
lattice = perm.lattice.d2_n8
saw_chain, weight = perm.mc_saw_rosenbluth_sampling(monomers=monomers, tries=mc_tries, lattice=lattice)
if saw_chain.monomers == 0:
    print("MC did not found a self avoiding walk.")
else:
    print(saw_chain)
    print(weight)

plot.plot_chain_2D(saw_chain)
```

    chain.monomers= 16
    0 0 0
    -1 0 0
    -1 -1 0
    -1 -2 0
    -2 -2 0
    -3 -1 0
    -2 0 0
    -2 1 0
    -1 2 0
    -1 1 0
    0 2 0
    -1 3 0
    0 4 0
    0 3 0
    1 3 0
    2 2 0
    
    313658956800.0



![png](scripts/notebooks/plot_chain_files/plot_chain_14_1.png)


### Draw many chains with rosenbluth samplings


```python
num_experiments = 10
monomers = 40
mc_tries = 100
dim = 2
lattice = perm.lattice.d2_n4
fig = plt.figure()
if dim == 2:
    ax = fig.add_subplot(111)
else:
    ax = fig.add_subplot(111, projection='3d')

saw_chains = []
weights = []
for _ in range(num_experiments):
    saw_chain, weight = perm.mc_saw_rosenbluth_sampling(monomers=monomers, tries=mc_tries, lattice=lattice)
    saw_chains.append(saw_chain)
    print(weight)
    weights.append(weight)
    if dim == 2:
        plot.plot_chain_2D(saw_chain, 0, 1, ax)
    else:
        plot.plot_chain(saw_chain, ax)
hist_weights, bins_weights = np.histogram(weights)
distro_weights = scipy.stats.rv_histogram([hist_weights, bins_weights])
fig = plt.figure()
ay = fig.add_subplot(111)
bins_weights_centers = (bins_weights[:-1] + bins_weights[1:])/2
ay.plot(bins_weights_centers, distro_weights.pdf(bins_weights_centers))
ay.set_xlabel('weights')
ay.set_ylabel('pdf')

```

    2.3425835473880064e+16
    2.3425835473880064e+16
    5.270812981623014e+16
    3.557798762595535e+17
    6169767367606272.0
    9254651051409408.0
    5205741216417792.0
    3.1624877889738086e+17
    1.0411482432835584e+16
    1.5617223649253376e+16





    Text(0, 0.5, 'pdf')




![png](scripts/notebooks/plot_chain_files/plot_chain_16_2.png)



![png](scripts/notebooks/plot_chain_files/plot_chain_16_3.png)


### PERM


```python
# perm uses recursion
num_experiments = 5
monomers = 400
mc_tries = 100000
# dim = 2
# lattice = perm.lattice.d2_n4
dim = 3
lattice = perm.lattice.d3_n26
k_bending = 5
energy_grow_func = perm.energy_grow_zero
# energy_grow_func = perm.energy_grow_bending(k=k_bending)
cube_half_side = 100
boundary_func = perm.boundary.is_inside_hyper_cube(half_side=cube_half_side)
half_side_x = 100
half_side_y = 100
half_side_z = 100
boundary_func = perm.boundary.is_inside_hyper_rectangle(half_side_x=half_side_x, half_side_y=half_side_y, half_side_z=half_side_z)
boundary_func = perm.boundary.is_inside_always
saw_chain, weight = perm.mc_saw_perm(monomers=monomers, tries=mc_tries,
                                     lattice=lattice, energy_grow_func=energy_grow_func,
                                     is_inside_boundary_func=boundary_func
                                     )
print("chain weight: ", weight)
print("ratio contour_length/ete_distance: ", saw_chain.contour_length()/saw_chain.ete_distance())
fig = plt.figure()
if dim == 2:
    ax = fig.add_subplot(111)
    plot.plot_chain_2D(saw_chain, 0, 1, ax)
else:
    ax = fig.add_subplot(111, projection='3d')
    plot.plot_chain(saw_chain, ax)

saw_chains = []
weights = []
for _ in range(num_experiments):
    saw_chain, weight = perm.mc_saw_perm(monomers=monomers, tries=mc_tries,
                                         lattice=lattice,
                                         energy_grow_func=energy_grow_func,
                                         is_inside_boundary_func=boundary_func)
    saw_chains.append(saw_chain)
    weights.append(weight)
    if dim == 2:
        plot.plot_chain_2D(saw_chain, 0, 1, ax)
    else:
        plot.plot_chain(saw_chain, ax)
```

    chain weight:  inf
    ratio contour_length/ete_distance:  21.39149283887122



![png](scripts/notebooks/plot_chain_files/plot_chain_18_1.png)



```python

# Estimate the partition function for a chain of n monomers
# <Z_m>_N = 1/N * sum(weight)
# The low and high thresholds, t_n and T_n are then chosen to be
# t_n = c*<Z_m>_N
# T_n = C*<Z_m>_N
# Where c, and C are fixed so C/c ~= 10

# First run regular Rosenbluth, with t_n = [0,...,0] and T_n = [inf, ..., inf]
# to compute a mean weight for the chain (at each length): Z_monomers
num_experiments = 10
mc_tries = 100
dim = 3
lattice = perm.lattice.d3_n26
k_bending = 3
energy_grow_func = perm.energy_grow_bending(k=k_bending)
boundary_func = perm.boundary.is_inside_always

max_monomers = 100
# The first weights are fixed and known
Z_monomers = [1.0, ]
for m in range(2, max_monomers):
    weights = []
    for _ in range(num_experiments):
        saw_chain, weight = perm.mc_saw_perm(monomers=m, tries=mc_tries, lattice=lattice,
                                             energy_grow_func=energy_grow_func,
                                             is_inside_boundary_func=boundary_func)

        weights.append(weight)
    Z_monomers.append(sum(weights) / num_experiments)

# Perform PERM with thresholds
weight_high_constant = 5
weight_low_constant = 0.5
weight_threshold_high = [weight_high_constant * i for i in Z_monomers]
weight_threshold_low = [weight_low_constant * i for i in Z_monomers]
saw_chains = []
weights = []
# for _ in range(num_experiments):
for _ in range(10):
    saw_chain, weight = perm.mc_saw_perm(monomers=max_monomers, tries=mc_tries,
                                         lattice=lattice,
                                         energy_grow_func=energy_grow_func,
                                         is_inside_boundary_func=boundary_func)
    saw_chains.append(saw_chain)
    weights.append(weight)

fig = plt.figure()
if dim == 2:
    ax = fig.add_subplot(111)
else:
    ax = fig.add_subplot(111, projection='3d')

for saw_chain in saw_chains:
    if dim == 2:
        plot.plot_chain_2D(saw_chain, 0, 1, ax)
    else:
        plot.plot_chain(saw_chain, ax)
```


![png](scripts/notebooks/plot_chain_files/plot_chain_19_0.png)

