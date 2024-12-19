# THEORY -- MAXWELL-BLOCH EQUATIONS FOR A MULTI-LEVEL RYDBERG ATOM

## Rotating frame and associated approximations

We generalise here the rotating frame method to a system of $N$ atomic levels having energies $\omega_i$ as dressed by fields described in README. Let the Rabi frequency of the $k\text{th}$ field coupling $|i\rangle$ and $|j\rangle$ be $\Omega_{ij,k}$. Setting $\hbar\rightarrow 1$ the total Hamiltonian is
$$
H = \sum_{i=1}^N \omega_i |i\rangle \langle i| + \sum_k \sum_{i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i \omega_k t} |i\rangle\langle j| + \text{h.c.}\right).
$$
Herein subscripts $i,j$ refer to atomic energy levels and a subscript $k$ refers to a laser field angular frequency.

We construct a rotating frame unitary transformation which strategically factors off driving laser angular frequencies from the various atomic transitions. Let
$$
|\psi(t)\rangle \rightarrow |\psi^\prime(t)\rangle = U(t) |\psi(t)\rangle,
$$
where
$$
U(t) = \sum_{i=1}^N e^{i \tilde{\omega}_i t} | i \rangle\langle i |
$$
for $\tilde{\omega}_i$ a set of rotating frame angular frequencies, which we set to
$$
\begin{align}
\tilde{\omega}_g &= \omega_g \\
\tilde{\omega}_e &= \omega_g + \omega_p \\
\tilde{\omega}_r &= \omega_e + \omega_c, r \in R
\end{align}
$$
where $R$ is the set of Rydberg states. The Hamiltonian is transformed to
$$
\begin{align}
H^\prime &= U H U^\dagger - i U \frac{dU^\dagger}{dt} \\
&= -\sum_i \delta_i |i\rangle\langle i| + \sum_{k,i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i(\omega_k - \tilde{\Delta}_{ij}) t} |i\rangle\langle j| + \text{h.c.}\right)
\end{align}
$$
where $\tilde{\Delta}_{ij} = \tilde{\omega}_i - \tilde{\omega}_j$.

We now drop fast oscillating terms; i.e., terms for which $\omega_k$ is far from $\tilde{\Delta}_{ij}$ such that the laser field is effectively decoupled from the transition. Physically, this requires that the laser bandwidth, the Rabi frequency, the Doppler broadening, and the natural linewidth together do not imply potential overlap. A typical laser bandwidth is on the order of 10 MHz, our Rabi frequencies are on the order of 10s of MHz, a Doppler broadening at room temperature is near 200 MHz at FWHM, and a natural linewidth on the order of 10 MHz.

Upon profiling the transitions ('transition_profiling.py'), we find that nearly all transitions deviate from the target 6s1/2--6p3/2 or 6p3/2--34d5/2 transitions by amounts on the order of THz, and thus do not couple at all to the laser fields; one single exception is the 6p3/2--34d3/2 transition, which differs from the 6p3/2--34d5/2 transition by only 2 GHz, putting it very slightly near the driving field upon accounting for Doppler broadening. Thus, we essentially need only retain Rabi terms on the three transitions
+ 6s1/2--6p3/2 (by design of the pump laser)
+ 6p3/2--34d5/2 (by design of the control laser)
+ 6p3/2--34d3/2 (incidental coupling to the control laser)

Let us denote these states of particular interest as
+ 6s1/2 : $|g\rangle$
+ 6p3/2 : $|e\rangle$
+ 34d5/2 : $|r1\rangle$
+ 34d3/2 : $|r2\rangle$

CONTINUE HERE:
Note to self: should introduce rotating frame factors on every single state; if not coupled to the laser field, then simply factor off their energies. If coupled to the laser field, proceed according to standard method (mostly completed here, needs some proofreading).


such that
$$
H^\prime = -\delta_p |e\rangle\langle e| + \frac{\Omega_{ge, p}}{2} \sigma^x_{ge} + \sum_{r \in R} \left(-\delta_{c,r} |r\rangle\langle r| + \frac{\Omega_{er, c}}{2} \sigma^x_{er}\right),
$$
with $\sigma^x_{ab}=|a\rangle\langle b| + |b\rangle\langle a|$.

## Dissipative terms and the master equation

We now move to a density operator picture in order to introduce dissipative and dephasing processes. We will introduce Doppler broadening in a later section. Let $\rho$ denote the density operator of the atomic system, and let $\mathcal{D}\left[L_{ab}\right]$ denote the dissipative superoperator describing decay from $|b\rangle$ to $|a\rangle$ with jump operator $L_{ab} = |a\rangle\langle b|$. $\mathcal{D}\left[L_{ab}\right]$ acts on $\rho$ according to
$$
\mathcal{D}\left[L_{ab}\right] \rho = L_{ab} \rho L^\dagger_{ab} - \frac{1}{2} \left(L^\dagger_{ab} L_{ab} \rho + \rho L^\dagger_{ab} L_{ab}\right).
$$
The association of all $\mathcal{D}[L_{ab}]$ with rates $\Gamma_{ab}$ defines the Lindblad superoperator
$$
\mathcal{L} \left[\rho\right] = \sum_{ab} \Gamma_{ab} \mathcal{D}\left[L_{ab}\right] \rho
$$
with which we assemble the master equation
$$
\frac{d \rho}{dt} = - i \left[H^\prime, \rho \right] + \mathcal{L}\left[\rho\right].
$$

Projecting the master equation onto the atomic basis, it has representation
$$
\begin{align}
\dot{\rho}_{gg} &= i \frac{\Omega_{ge,p}}{2} \left(\rho_{ge} - \rho_{eg}\right) + \Gamma_{ge} \rho_{ee} + \sum_r \Gamma_{gr} \rho_{rr} \\
\dot{\rho}_{eg} &= i \delta_p \rho_{eg} + i \frac{\Omega_{ge,p}}{2} \left(\rho_{ee} - \rho_{gg}\right) - \frac{\Gamma_{ge}}{2} \rho_{eg} \\
\dot{\rho}_{rg} &= i\delta_{c,r} \rho_{rg} - \frac{\Gamma_{gr}}{2} \rho_{rg}\\
\dot{\rho}_{ee} &= i \frac{\Omega_{ge,p}}{2} \left(\rho_{eg} - \rho_{ge}\right) + i \sum_r \frac{\Omega_{er,c}}{2}\left(\rho_{er}-\rho_{re}\right) - \Gamma_{ge} \rho_{ee} + \sum_r \Gamma_{er} \rho_{rr} \\
\dot{\rho}_{re} &= i \left(\delta_{c,r}-\delta_p\right) \rho_{re} + i \frac{\Omega_{er,c}}{2} \left(\rho_{rr}-\rho_{ee}\right) - \frac{\Gamma_{er}}{2} \rho_{re} \\
\dot{\rho}_{rr} &= i \frac{\Omega_{er,c}}{2}\left(\rho_{re} - \rho_{er} \right) - \left(\Gamma_{er} + \Gamma_{gr}\right) \rho_{er}
\end{align}
$$
