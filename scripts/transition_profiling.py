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
import time

# Append parent to path for resolving imports in adjacent folders
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Imports from adjacent folders
from fns.fns_rydberg_mb import *
from classes.classes_rydberg_mb import *

# Sim size
nmin = 6  # Minimum n
nmax = 35  # Maximum n
lmin = 0  # Minimum l
lmax = 5  # Maximum l

# Constants and conversions
c_light = 2.99792e8             # m / s
h_planck = 6.62607e-34          # SI, J / Hz
hbar_planck = h_planck / (2. * np.pi)   # SI, J / (rad/s)
e_coulomb = 1.60218e-19         # C
a0_metres = 5.29177e-11         # m
eps0_si = 8.854187817e-12       # C^2 / (J m)
k_boltz = 1.3806e-23            # Boltzmann constant (J/K)
t_kelvin = 298                  # Temperature in Kelvin (25 deg C)
m_kg = 2.21e-25                 # Mass of caesium atom in kg

joules_to_mhz = 1e-6 / h_planck
ev_to_mhz = e_coulomb * joules_to_mhz
ea0_to_mhz_v_m = e_coulomb * a0_metres * joules_to_mhz
v_fwhm = 1.177 * np.sqrt(2. * k_boltz * t_kelvin / m_kg) 

# Load parameters for Caesium
atom = Caesium()

""" ASSEMBLE THE ENERGY LEVELS
Pull energy levels from the ARC database
 and append to a dictionary with key elts
 corresponding to state strings

"""

levels_list = []
index_list = []
for this_state in AtomicNLJMIterator(nmin, nmax, lmax, iter_mj=False):
    this_energy = atom.getEnergy(*this_state) * ev_to_mhz
    levels_list.append({'energy':this_energy})
    index_list.append(this_state)
    
levels_df = pd.DataFrame(levels_list, index=index_list)

print(levels_df)

""" ASSEMBLE TRANSITIONS

Assemble a dataframe of the transitions,
 indexed by 'spec_name_a-spec_name_b',
 including:
    'energy' : energy in MHz
    'dip' : dipole elt in MHz/(V/m)

"""
trans_list = []
index_list = []
for state_a in AtomicNLJMIterator(nmin,nmax,lmax):
    for state_b in AtomicNLJMIterator(nmin,nmax,lmax):
        if not is_lin_dip_permitted(state_a, state_b) : continue     # Store only dipole permitted transitions
        if state_a[3] != .5 : continue      # Wlog choose only mj=0.5 states
        if state_b[1] <= state_a[1] : continue      # l1 < l2 by convention for transitions
        dip_mhz_v_m = atom.getDipoleMatrixElement(*state_a, *state_b, 0) * ea0_to_mhz_v_m
        delta_mhz = levels_df.at[state_b[:3], 'energy'] - levels_df.at[state_a[:3], 'energy']
        omega_si = delta_mhz * 1e6 * 2. * np.pi
        dip_si = dip_mhz_v_m / joules_to_mhz
        gamma_us = 1e-6 * (
            omega_si ** 3 * dip_si ** 2
            / (3. * np.pi * eps0_si * hbar_planck * c_light ** 3))
        lwid_mhz = gamma_us / (2. * np.pi)
        doppler_mhz = delta_mhz * v_fwhm / c_light
        trans_tuple = (*state_a[:3], *state_b[:3])
        trans_list.append({
            'delta' : delta_mhz,
            'dip' : dip_mhz_v_m,
            'gamma' : abs(gamma_us),
            'lwid' : abs(lwid_mhz),
            'doppler' : doppler_mhz})
        index_list.append(trans_tuple)
trans_df = pd.DataFrame(trans_list, index=index_list)

print('The transitions are')
print(trans_df)

print('The transition deviations from 6s1/2--6p3/2 are (ranked):')
print(compare_transitions((6,0,0.5,6,1,1.5), trans_df))
print('-')

print('The transition deviations from 6p3/2--34d5/2 are (ranked):')
print(compare_transitions((6,1,1.5,34,2,2.5), trans_df))
print('-')

