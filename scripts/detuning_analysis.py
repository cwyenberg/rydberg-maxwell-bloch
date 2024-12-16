""" detuning_analysis.py

Analyze the detuning of the Caesium levels from the driving fields
in order to quantify the degree of coupling to unintended transitions.

The QVIL pump field acts on the 6S1/2 to 6P3/2 transition near 852 nm,
and the control field on the 6P3/2 to 34D5/2 transition at 509 nm.

"""

import numpy as np
from matplotlib import pyplot as plt
from arc import *
from itertools import product

# Constants
c_light = 2.99792e8
h_planck = 6.62607e-34
e_coulomb = 1.60218e-19
ev_to_mhz = 1e-6 * e_coulomb / h_planck

# Load parameters for Caesium
atom = Caesium()

nmin = 6  # Minimum n
nmax = 34  # Maximum n
lmin = 0  # Minimum l
lmax = 3  # Maximum l

# Quick function to build a state string
#  from quantum numbers
def build_state_str(n:int, l:int, j:float):
    orbs = ['s','p','d','f','g','h']
    return str(n) + orbs[l] + str(j)

# Get the energy levels:
#  Pull energy levels from the ARC database
#  and append to a dictionary with key elts
#  corresponding to state strings
lvls = {}
for n in range(nmin, nmax+1):
    for l in range(lmin, lmax+1):
        jmin = max(l-.5, .5)
        jmax = l+.5
        for j in np.arange(jmin, jmax+1):
            this_energy = atom.getEnergy(n,l,j)
            this_state_str = build_state_str(n,l,j)
            lvls[this_state_str] = this_energy

# Get energy deltas:
lvl_deltas = {}
#  Enumerate all energy differences between non-forbidden
#  transitions
nrange = range(nmin, nmax+1)
lrange = range(lmin, lmax+1)
spin_comb_range = [-.5,.5]
delta_j_range = [0., 1.]
# The a state ranges over all quantum numbers, but one less
#  than the max l: the state a is considered the lower l state.
# We enumerate over j = l+/-1/2 corresponding to spin combination
for (n_a, l_a, jos_a) in product(
    range(nmin, nmax+1), range(lmin, lmax), spin_comb_range):
    
    j_a = l_a + jos_a
    if j_a < 0 : continue       # Skip if j_a < 0

    state_a_str = build_state_str(n_a,l_a,j_a)

    l_b = l_a + 1
    for (n_b, delta_j) in product(
        range(nmin, nmax+1), delta_j_range):

        j_b = j_a + delta_j
        if j_b < l_b-.5 or j_b > l_b+.5 : continue   # Skip if j_b OoR

        state_b_str = build_state_str(n_b,l_b,j_b)

        this_key_str = state_a_str + '-' + state_b_str

        lvl_deltas[this_key_str] = lvls[state_b_str] - lvls[state_a_str]

pump_line_mhz = lvl_deltas['6s0.5-6p1.5'] * ev_to_mhz
pump_line_nm = (c_light / pump_line_mhz) * 1e3
print(pump_line_nm)

ctrl_line_mhz = lvl_deltas['6p1.5-34d2.5'] * ev_to_mhz
ctrl_line_nm = (c_light / ctrl_line_mhz) * 1e3
print(ctrl_line_nm)


