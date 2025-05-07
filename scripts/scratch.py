import numpy as np
import matplotlib.pyplot as plt

# Parameters
tmax = 1000.
nt = 10000
time = np.linspace(0., tmax, nt, endpoint=False)
frequency = 1.  # Hz
lw = .1
num_oscillations = 100

# Generate superposition of oscillating exponentials
superposition = np.zeros_like(time, dtype=complex)
for _ in range(num_oscillations):
    phase = np.random.uniform(0., 2. * np.pi)
    frequency = np.random.normal(frequency, lw)
    superposition += np.exp(1j * (2 * np.pi * frequency * time + phase))

# Plot the real part of the superposition
plt.plot(time, superposition.real)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Superposition of Oscillating Exponentials')
plt.show()
