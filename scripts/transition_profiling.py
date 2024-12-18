""" Transition profiling

Analyze the detuning of the Caesium levels from the driving fields
in order to quantify the degree of coupling to unintended transitions.

The QVIL pump field acts on the 6S1/2 to 6P3/2 transition near 852 nm,
and the control field on the 6P3/2 to 34D5/2 transition at 509 nm.

ON PERMITTED TRANSITIONS
On the topic of permitted transitions, there are multiple factors to consider:
1. Conventional conservation laws: Delta j = 0, +/-1; Delta mj = 0, +/- 1
2. Specific symmetries of the dipole transition (which dominates): the wavefunction
must change parity; i.e., Delta l = +/- 1; and the interaction does not involve
the magnetic moment; i.e., mj will not change.
Note that Delta j = 0 can occur; this does not contradict conservation of total
angular momentum due to the following subtlety: although the photon carries away
1 hbar, j is related to the magnitude of the total angular momentum, and not the
vector quantity; so it may remain unchanged (think, classically, of a vector
swapping direction).

Practically, the only dipole transitions we witness obey:

Delta l = +/-1
Delta j = 0, +/-1
Delta mj = 0

UNIT CONVENTIONS
Energy in MHz
Dipole in MHz/(V/m)
Gamma in us^-1 (note not referred to as MHz, since appears as coeff in diff eq)

"""

# Generic imports
import numpy as np
from matplotlib import pyplot as plt
from arc import *
from itertools import product
import pandas as pd

# Append parent to path for resolving imports in adjacent folders
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Imports from adjacent folders
from fns.fns_rydberg_mb import *
from classes.classes_rydberg_mb import *

# Sim size
nmin = 6  # Minimum n
nmax = 7  # Maximum n
lmin = 0  # Minimum l
lmax = 1  # Maximum l

# Constants and conversions
c_light = 2.99792e8             # m / s
h_planck = 6.62607e-34          # SI, J / Hz
hbar_planck = h_planck / (2. * np.pi)   # SI, J / (rad/s)
e_coulomb = 1.60218e-19         # C
a0_metres = 5.29177e-11         # m
eps0_si = 8.854187817e-12       # C^2 / (J m)

joules_to_mhz = 1e-6 / h_planck
ev_to_mhz = e_coulomb * joules_to_mhz
ea0_to_mhz_v_m = e_coulomb * a0_metres * joules_to_mhz

# Load parameters for Caesium
atom = Caesium()

""" ASSEMBLE THE ENERGY LEVELS
Pull energy levels from the ARC database
 and append to a dictionary with key elts
 corresponding to state strings
"""

levels_df = pd.DataFrame()
for state_tuple in AtomicNLJMIterator(nmin, nmax, lmax, iter_mj=True):
    this_energy = atom.getEnergy(*state_tuple[:3]) * ev_to_mhz
    this_entry = pd.DataFrame({
        'energy' : this_energy},
        index = [state_tuple])
    levels_df = pd.concat([levels_df, this_entry])

print(levels_df)

""" ASSEMBLE TRANSITIONS

Assemble a dataframe of the transitions,
 indexed by 'spec_name_a-spec_name_b',
 including:
    'energy' : energy in MHz
    'dip' : dipole elt in MHz/(V/m)

"""
trans_df = pd.DataFrame()
for state_a in AtomicNLJMIterator(nmin,nmax,lmax):
    for state_b in AtomicNLJMIterator(nmin,nmax,lmax):
        if not is_lin_dip_permitted(state_a, state_b): continue     # Store only dipole permitted transitions
        if state_b[1] <= state_a[1] : continue      # l1 < l2 by convention for transitions
        trans_tuple = (*state_a, *state_b)
        dip = atom.getDipoleMatrixElement(*state_a, *state_b, 0) * ea0_to_mhz_v_m
        delta = levels_df.at[state_a,'energy'] - levels_df.at[state_b,'energy']
        omega_si = delta * 1e6 * 2. * np.pi
        dip_si = dip / joules_to_mhz
        gamma_us = 1e-6 * (
            omega_si ** 3 * dip_si ** 2
            / (3. * np.pi * eps0_si * hbar_planck * c_light ** 3))
        this_entry = pd.DataFrame({
            'delta' : delta,
            'dip' : dip,
            'gamma' : gamma_us},
            index = [trans_tuple])
        trans_df = pd.concat([trans_df, this_entry])

print(trans_df)

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
        this_trans_str = state_a_str + '-' + state_b_str
        trans_delta = (
            levels_df.loc[state_b_str,'energy']
            - levels_df.loc[state_a_str,'energy'])
        this_entry = pd.DataFrame({
            'delta' : trans_delta},
            index = [this_trans_str])
        trans_df = pd.concat([trans_df, this_entry])


# DIPOLE ELEMENTS
for trans_name, row_value in trans_df.iterrows():
    # print(trans_name)
    # print(row_value)
    left_state, right_state = trans_name.split('-')
    print(left_state)
    print(right_state)
    atom.getDipoleMatrixElement()



# pump_line_mhz = lvl_deltas_mhz['6s0.5-6p1.5']
# pump_line_nm = (c_light / pump_line_mhz) * 1e3
# print('The pump drives at ' + str(pump_line_mhz) + ' MHz')
# # print(pump_line_nm)

# # ctrl_line_mhz = lvl_deltas_mhz['6p1.5-34d2.5']
# # ctrl_line_nm = (c_light / ctrl_line_mhz) * 1e3
# # print('The control drives at ' + str(ctrl_line_mhz) + ' MHz')
# # print(ctrl_line_nm)

# def find_nearest_trans(target_trans_key:str, lvl_deltas_mhz:dict,num=2):
#     """ Find the nearest transition to a target transition """
#     target_trans_mhz = abs(lvl_deltas_mhz[target_trans_key])
#     print('')
#     print('The target transition is ' + target_trans_key)
#     print(' of energy difference ' + str(target_trans_mhz) + ' MHz.')
#     print('')
#     print('The nearest ' + str(num) + ' transitions have the following sorted distances:')

#     other_trans_df = pd.DataFrame()

#     for key, value in lvl_deltas_mhz.items():
#         this_diff_mhz = abs(value - target_trans_mhz)
#         this_entry = pd.DataFrame({
#             'trans energy': value,
#             'difference' : this_diff_mhz},
#             index = [key])
#         other_trans_df = pd.concat([other_trans_df,this_entry])

#     other_trans_df = other_trans_df.sort_values(by='difference')

#     print(other_trans_df[:num])

# find_nearest_trans('6s0.5-6p1.5', lvl_deltas_mhz,num=5)
# find_nearest_trans('6p1.5-34d2.5', lvl_deltas_mhz,num=5)
