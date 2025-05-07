""" Sets up a Rydberg atom """

from arc import Cesium, PairStateInteractions

# Define the Cesium atom
cs = Cesium()

# Define the principal quantum number and state
n = 34
l = 2  # D-state (l = 2)
j = 2.5  # Total angular momentum J = 5/2
mj1 = .5
mj2 = .5

state34d5o2 = (n, l, j)

# Compute the van der Waals coefficient (C6)
C6_coefficient = cs.getC6term(*state34d5o2, *state34d5o2, *state34d5o2)

# Print the result
print(f"C6 coefficient for Cs(34D5/2) pair-state: {C6_coefficient:.2e} GHz (µm)^6")
