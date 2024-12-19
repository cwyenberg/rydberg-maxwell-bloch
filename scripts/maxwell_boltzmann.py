import numpy as np
import matplotlib.pyplot as plt

""" MAXWELL-BOLTZMANN DISTRIBUTION

Investigate the MB distribution for Caesium atoms
and its effect upon linewidths

"""

# Constants
k_B = 1.3806e-23  # Boltzmann constant (J/K)
T = 298  # Temperature in Kelvin (25 deg C)
m = 2.21e-25  # Mass of caesium atom in kg

# Maxwell-Boltzmann distribution function
def maxwell_boltzmann(v, T, m):
    prefactor = (m / (2 * np.pi * k_B * T))**(3/2)
    return prefactor * v**2 * np.exp(-m * v**2 / (2 * k_B * T))

# Speed range
v = np.linspace(0, 1000, 500)  # Speed range in m/s

# Compute the Maxwell-Boltzmann distribution for caesium
f_v = maxwell_boltzmann(v, T, m)

v_fwhm = 1.177 * np.sqrt(2. * k_B * T / m) 
print('The FWHM at this temperature is ' + str(v_fwhm) + ' m/s.')
print('plotting...')

# Plot the distribution
plt.plot(v, f_v, label='Maxwell-Boltzmann Distribution', color='blue')
plt.title('Maxwell-Boltzmann Distribution for Caesium at Room Temperature')
plt.xlabel('Speed (m/s)')
plt.ylabel('Probability Density')
plt.grid(True)
plt.legend()
plt.show()
