# Maxwell-Bloch Model of a Caesium vapour atomic sensor

We consider a thermal ensemble of Caesium atoms enclosed within a vapour cell and dressed with multiple laser fields. The laser fields are gaussian beams of geometries detailed where they appear within this project.

We treat the system semi-classically, where the laser fields modify the system Hamiltonian living in the Hilbert space of atomic states but not referencing photonic number states. This approximation is valid in the strong-field limit, which is typically the case for atomic sensors.

The laser fields possess such small bandwidth (much less than any level separation) that we neglect their spectral extent; however, the thermal distribution broadens the spectrum such that multiple Rydberg states must be included (the precise number depending upon the temperature). The lower atomic levels are, however, sufficiently separated such that intermediate levels between those on resonance with the lasers need not be modelled; however, each lower level is also broadened (which has a strong effect upon laser coupling).

The states of the system are described by
+ A thermal distribution $P_v(T)$, corresponding to velocity classes and corresponding density operators $\rho_v$.
+ A spatial volume corresponding to classical atomic coordinates, necessitated only by the non-uniform Gaussian beam
+ Atomic levels chosen by consideration of the dressing fields and Doppler broadening extent
+ An intermediate "bath state" necessitated by spontaneous emission

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

