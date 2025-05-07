from arc import Cesium, PairStateInteractions

cs = Cesium()

# Dipole-Interaction Dispersion Coefficient:60S1/2
# ================================================
# Evaluation of the Cs 30D_5/2 C6 coefficient using perturbation theory (Theta=0,phi=0)
n0 = 40
l0 = 2
j0 = 2.5
mj0 = 2.5
# Target State
theta = 0
# Polar Angle [0-pi]
phi = 0
# Azimuthal Angle [0-2pi]
dn = 5
# Range of n to consider (n0-dn:n0+dn)
deltaMax = 25e9  # Max pair-state energy difference [Hz]

# Set target-state and extract value
calculation = PairStateInteractions(
    cs, n0, l0, j0, n0, l0, j0, mj0, mj0
)
C6 = calculation.getC6perturbatively(theta, phi, dn, deltaMax)
print("C6 [%s] = %.2f GHz (mum)^6" % ((n0, l0, j0), C6))
