# Maxwell-Bloch Model of a Caesium vapour atomic sensor

We consider a thermal ensemble of Caesium atoms enclosed within a vapour cell and dressed with multiple laser fields. The laser fields are gaussian beams of geometries detailed where they appear within this project.

We treat the system semi-classically, where the laser fields modify the system Hamiltonian living in the Hilbert space of atomic states but not referencing photonic number states. This approximation is valid in the strong-field limit which is typically the case for atomic sensors.

The laser fields possess such small bandwidth (much less than any level separation) that we neglect their spectral extent; however, Doppler broadening and spontaneous emission processes require that multiple atomic states be included. The atomic level spacing and broadening are, as analyzed elsewhere, of sufficient structures such that the lasers couple non-negligibly only to the engineered transitions; however, the broadening of the relevant levels has a strong effect upon the coupling strength.

The state of the system is described by
+ A thermal distribution $P_v(T)$, corresponding to velocity classes and corresponding density operators $\rho_v$.
+ Atomic levels chosen by consideration of the dressing fields and Doppler broadening extent

## Level structure

The QVIL Caesium vapour cell operates between the following atomic states:
+ $R$ : Additional rydberg states $|r^\prime\rangle, \omega_{r^\prime}$
+ $34D_{5/2}$ : Rydberg state $|r\rangle, \omega_r$
+ $6P_{3/2}$ : Excited state $|e\rangle, \omega_e$
+ $6S_{1/2}$ : Ground state $|g\rangle, \omega_g$

The pump field of Rabi strength $\Omega_p$ has angular frequency $\omega_p = \omega_{ge} + \delta_{p,eg}$ near the $g,e$ transition. The control field $\Omega_c$ has $\omega_c = \omega_{er} + \delta_{c,er}$ near the $e,r$ transition. QVIL operates its pump laser near an 852 nm $g,e$ transition and its control laser near a 511 nm $e,r$ transition.

## Steady-state considerations

In steady-state conditions, $\frac{d\rho(v)}{dt} = 0$, we solve the equations for each velocity class assuming steady-state populations and coherences. The steady-state density matrix for the Doppler-broadened system is then:

$$
\rho_{\text{total}} = \int_{-\infty}^{\infty} \rho_{\text{ss}}(v) f(v) \, dv.
$$

In the case of a two-level system, the Doppler-broadened spectrum is obtained by calculating the imaginary part of $\rho_{ab}$ [this remark to be generalized to the current many-level discussion].

