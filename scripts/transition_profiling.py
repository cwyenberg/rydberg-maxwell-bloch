""" Transition profiling

Analyze the detuning of the Caesium levels from the driving fields
in order to quantify the degree of coupling to unintended transitions.

The QVIL pump field acts on the 6S1/2 to 6P3/2 transition near 852 nm,
and the control field on the 6P3/2 to 34D5/2 transition at 509 nm.

ON PERMITTED TRANSITIONS
On the topic of permitted transitions, there are multiple factors to consider:
1. Conventional conservation laws: Delta j = 0, +/-1; Delta mj = 0, +/- 1
2. Specific symmetries of the dipole transition (which dominates): the wavefunction
must change parity; i.e., Delta l = +/- 1.
Note that Delta j = 0 can occur; this does not contradict conservation of total
angular momentum due to the following subtlety: although the photon carries away
1 hbar, j is related to the magnitude of the total angular momentum, and not the
vector quantity; so it may remain unchanged (think, classically, of a vector
swapping direction).

Practically, the only dipole transitions we witness obey:

Delta l = +/-1
Delta j = 0, +/-1
Delta mj = 0, +-1 (0 for linear, +/- 1 for left- or right-circular, unsure of respectively(?))

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

pd.set_option('display.max_rows', None)

# Sim size
nmin = 6  # Minimum n
nmax = 70  # Maximum n
lmin = 0  # Minimum l
lmax = 1  # Maximum l

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

""" DIAGNOSTICS """

state_a = (6, 1, 1.5, .5)
state_b = (6, 2, 1.5, -.5)
print('This transition has dip ')
print(atom.getDipoleMatrixElement(*state_a, *state_b, -1) * ea0_to_mhz_v_m)



print('Assembling energy levels...')
levels_list = []
index_list = []
for this_state in AtomicNLJMIterator(nmin, nmax, lmax, iter_mj=False):
    this_energy = atom.getEnergy(*this_state) * ev_to_mhz
    levels_list.append({'energy':this_energy})
    index_list.append(this_state)
    
levels_df = pd.DataFrame(levels_list, index=index_list)


""" ASSEMBLE TRANSITIONS

Assemble a dataframe of the transitions,
 indexed by 'spec_name_a-spec_name_b',
 including:
    'energy' : energy in MHz
    'dip' : dipole elt in MHz/(V/m)

"""

print('Assembling transitions...')
trans_list = []
index_list = []
for state_a in AtomicNLJMIterator(nmin,nmax,lmax):
    for state_b in AtomicNLJMIterator(nmin,nmax,lmax):
        dip_valid, dip_type = dip_trans_type(state_a, state_b)
        if not dip_valid : continue     # Store only dipole permitted transitions
        if state_b[1] <= state_a[1] : continue      # l1 < l2 to avoid double counting
        dip_mhz_v_m = atom.getDipoleMatrixElement(*state_a, *state_b, dip_type) * ea0_to_mhz_v_m
        delta_mhz = levels_df.at[state_b[:3], 'energy'] - levels_df.at[state_a[:3], 'energy']
        omega_si = delta_mhz * 1e6 * 2. * np.pi
        dip_si = dip_mhz_v_m / joules_to_mhz
        gamma_us = 1e-6 * abs(
            omega_si ** 3 * dip_si ** 2
            / (3. * np.pi * eps0_si * hbar_planck * c_light ** 3))
        lwid_mhz = gamma_us / (2. * np.pi)
        doppler_mhz = delta_mhz * v_fwhm / c_light
        trans_tuple = (*state_a, *state_b)
        trans_list.append({
            'delta' : delta_mhz,
            'dip' : dip_mhz_v_m,
            'dip-type' : dip_type,
            'gamma' : gamma_us,
            'lwid' : lwid_mhz,
            'doppler' : doppler_mhz})
        index_list.append(trans_tuple)
trans_df = pd.DataFrame(trans_list, index=index_list)

trans_lvls = np.arange(7,71)
trans_rates = np.empty(trans_lvls.shape[0], dtype=float)
trans_dips = np.empty(trans_lvls.shape[0], dtype=float)
trans_deltas = np.empty(trans_lvls.shape[0], dtype=float)
for id, top_lvl in enumerate(trans_lvls):
    trans_tuple = (6, 0, 0.5, 0.5, top_lvl, 1, 1.5, 0.5)
    trans_rates[id] = trans_df.at[trans_tuple, 'gamma']
    trans_dips[id] = trans_df.at[trans_tuple, 'dip']
    trans_deltas[id] = trans_df.at[trans_tuple, 'delta']

plt.cla()
plt.loglog(trans_lvls, trans_rates)
plt.title('Transition rates to ground vs principal')
plt.show()

plt.cla()
plt.loglog(trans_lvls, np.abs(trans_dips))
plt.title('Transition dipoles to ground vs principal')
plt.show()

plt.cla()
plt.loglog(trans_lvls, trans_deltas)
plt.title('Transition deltas to ground vs principal')
plt.show()


# plt.cla()
# print('Plotting the spontaneous emission rates...')
# gamma_vals = trans_df['gamma'].to_numpy()
# plt.scatter(range(gamma_vals.shape[0]),gamma_vals)
# plt.show()


print('The 6p3/2 to 34d5/2 transition Doppler broadening fwhm is ')
print(trans_df.at[(6,1,1.5,34,2,2.5), 'doppler'])

print('The 6s1/2 to 6p3/2 wavelength is ')
print(frequency_to_wavelength(trans_df.at[(6,0,0.5,6,1,1.5), 'delta']))
print('The top five nearest deviations from 6s1/2--6p3/2 are (ranked):')
ref_delta = trans_df.at[(6,0,0.5,6,1,1.5), 'delta']
print(compare_transitions(ref_delta, trans_df).head())
print('-')

print('The 6p3/2 to 34d5/2 wavelength is ')
print(frequency_to_wavelength(trans_df.at[(6,1,1.5,34,2,2.5), 'delta']))
print('The top five nearest deviations from 6p3/2--34d5/2 are (ranked):')
ref_delta = trans_df.at[(6,1,1.5,34,2,2.5), 'delta']
print(compare_transitions(ref_delta, trans_df).head())
print('-')

print('The top five nearest transitions to 852 nm are (ranked):')
ref_delta = wavelength_to_frequency(852)
print(compare_transitions(ref_delta, trans_df).head())
print('-')

print('The top five nearest transitions to 509 nm are (ranked):')
ref_delta = wavelength_to_frequency(509)
print(compare_transitions(ref_delta, trans_df).head())
print('-')
