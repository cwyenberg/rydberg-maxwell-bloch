"""

Generation of a Caesium EIT plot for review paper

C Wyenberg
May 7, 2025

"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from arc import Caesium

"""
QVIL reference values from presentation:
Omega_p evaluates to 300 kHz
Omega_c evaluates to 10 MHz

The levels are
|6S1/2, F=4, mF=0> = |g>
|6P3/2, F=5, mF=0> = |e>
|34D5/2, F=5, mF=0> = |r>

The coupling laser is on the |e> -> |r> transition
with a transition wavelength of 509 nm,
corresponding to a frequency of 589 THz,
though this value is not needed in the rotating frame.

The probe laser is on the |g> -> |e> transition
with a transition wavelength of 852 nm,
corresponding to a frequency of 352 THz,
which is also not needed in the rotating frame.

The physical decay rates are computed by ARC.

"""

# User config entries:
Detune_p_max = 80.0            # Probe detuning plot limit (MHz)
Rabi_p_weak = 0.5              # Weaker probe Rabi frequency (MHz)
Rabi_p_strong = 1.0            # Stronger probe Rabi frequency (MHz)
Rabi_c = 10.0                  # Coupling Rabi frequency (MHz)
Detune_c = 0.0                 # Coupling laser detuning (MHz)
n_Detune = 1000                # Number of probe detuning points

# Construct a Caesium atom instance from ARC
cs = Caesium()

# Convert to angular frequencies
Omega_p_weak = 2. * np.pi * Rabi_p_weak
Omega_p_strong = 2. * np.pi * Rabi_p_strong
Omega_c = 2. * np.pi * Rabi_c
Delta_c = 2. * np.pi * Detune_c

# Effective lifetimes (excluding blackbody radiation) in us
e_lifetime = 1e6 * cs.getStateLifetime(6, 1, 1.5)
r_lifetime = 1e6 * cs.getStateLifetime(34, 2, 2.5)

# Obtain system parameters from ARC
Gamma_eg = 2. * np.pi / e_lifetime              # Decay rate Gamma value from |e> to |g> (rad/us)
Gamma_re = 2. * np.pi / r_lifetime              # Decay rate Gamma value from |r> to |e> (rad/us)

# Basis states: |g>, |e>, |r>
ket_g = basis(3, 0)
ket_e = basis(3, 1)
ket_r = basis(3, 2)

# Operators
sigma_ge = ket_g * ket_e.dag()
sigma_eg = ket_e * ket_g.dag()
sigma_er = ket_e * ket_r.dag()
sigma_re = ket_r * ket_e.dag()

# Collapse operators
C1 = np.sqrt(Gamma_eg) * sigma_ge
C2 = np.sqrt(Gamma_re) * sigma_er
collapse_ops = [C1, C2]

# Detuning range for probe
Delta_p_vals = 2. * np.pi * np.linspace(-Detune_p_max, Detune_p_max, n_Detune)   # (rad/us)

# Storage for probe absorption (Im[ρ_ge])
absorption_weak = []
absorption_strong = []

# Loop over probe detunings
for Delta_p in Delta_p_vals:

    # Hamiltonians in rotating frame
    H_weak = (
        .5 * Omega_p_weak * (sigma_eg + sigma_ge)
        + .5 * Omega_c * (sigma_re + sigma_er)
        - Delta_p * ket_e * ket_e.dag()
        - (Delta_p + Delta_c) * ket_r * ket_r.dag()
    )
    H_strong = (
        .5 * Omega_p_strong * (sigma_eg + sigma_ge)
        + .5 * Omega_c * (sigma_re + sigma_er)
        - Delta_p * ket_e * ket_e.dag()
        - (Delta_p + Delta_c) * ket_r * ket_r.dag()
    )

    # Solve for steady states
    rho_ss_weak = steadystate(H_weak, collapse_ops)
    rho_ss_strong = steadystate(H_strong, collapse_ops)

    # Absorptions ~ Im[ρ_ge]
    absorption_weak.append(np.abs(np.imag((rho_ss_weak * sigma_ge).tr())))
    absorption_strong.append(np.abs(np.imag((rho_ss_strong * sigma_ge).tr())))

# Plotting
plt.rcParams.update({'font.size': 10})
plt.figure(figsize=(6, 4))
plt.plot(Delta_p_vals / (2 * np.pi), absorption_weak, label='Probe Rabi frequency 0.5 MHz', color="#CF123D")
plt.plot(Delta_p_vals / (2 * np.pi), absorption_strong, label='Probe Rabi frequency 1.0 MHz', color="#0492D2")
plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
plt.xlabel("Probe Detuning (MHz)")
plt.ylabel("Probe Absorption (arb.)")
# plt.title("Electromagnetically Induced Transparency (EIT) in Caesium Vapour Cell")
# plt.grid(True)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('caesium_eit.png', dpi=300)
plt.savefig('caesium_eit.eps')
plt.show()
