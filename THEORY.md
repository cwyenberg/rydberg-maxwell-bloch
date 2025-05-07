# THEORY
# MAXWELL-BLOCH EQUATIONS FOR A MULTI-LEVEL RYDBERG ATOM

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
for $\tilde{\omega}_i$ a set of rotating frame angular frequencies to be defined shortly. The Hamiltonian is transformed to
$$
\begin{align}
H^\prime &= U H U^\dagger - i U \frac{dU^\dagger}{dt} \\
&= -\sum_i \delta_i |i\rangle\langle i| + \sum_{k,i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i(\omega_k - \tilde{\Delta}_{ij}) t} |i\rangle\langle j| + \text{h.c.}\right)
\end{align}
$$
where $\delta_i = \tilde{\omega}_i - \omega_i$ and $\tilde{\Delta}_{ij} = \tilde{\omega}_i - \tilde{\omega}_j$.

We now drop fast oscillating terms; i.e., terms for which $\omega_k$ is `far' from $|\omega_i - \omega_j|$ such that the laser field is effectively decoupled from the transition. More precisely, 'far' implies that the laser bandwidth, the Rabi frequency, the Doppler broadening, and the natural linewidth together do not present potential overlap. A typical laser bandwidth is on the order of 10s of MHz, our Rabi frequencies are on the order of 10s of MHz, a Doppler broadening at room temperature is near 200 MHz at FWHM, and a natural linewidth on the order of 10s of MHz.

Upon profiling the transitions ('transition_profiling.py'), we find that nearly all transitions deviate from the target 6s1/2--6p3/2 or 6p3/2--34d5/2 transitions by amounts on the order of THz, and therefore do not couple to the laser fields; one single exception is the 6p3/2--34d3/2 transition, which differs from the 6p3/2--34d5/2 transition by only 2 GHz, putting it very slightly near the driving field upon accounting for Doppler broadening. Thus, we essentially need only retain Rabi terms on the three transitions
+ 6s1/2--6p3/2: (by design of the pump laser),
+ 6p3/2--34d5/2 (by design of the control laser), and
+ 6p3/2--34d3/2 (incidental coupling to the control laser).

Let us denote these states of particular interest as
+ 6s1/2 : $|g\rangle$,
+ 6p3/2 : $|e\rangle$,
+ 34d5/2 : $|r_1\rangle$, and
+ 34d3/2 : $|r_2\rangle$.

We will herein refer to the above four states ($g,e,r_1,r_2$) as the primary states, and all other states as the peripheral states. Having identified the three non-trivial transitions ${g \leftrightarrow e}, {e \leftrightarrow r_1},$ and ${e\leftrightarrow r_2}$, we may prescribe convenient rotating frame factors 
$$
\begin{align}
\tilde{\omega}_g &= \omega_g, \\
\tilde{\omega}_e &= \omega_g + \omega_p, \\
\tilde{\omega}_{r_1/r_2} &= \omega_e + \omega_c, \text{ and} \\
\tilde{\omega}_a &= \omega_a \text{ ($a$ denotes all remaining states).}
\end{align}
$$
such that
$$
H^\prime = -\delta_p |e\rangle\langle e| - \delta_{r_1} |r_1\rangle\langle r_1| - \delta_{r_2} |r_2\rangle\langle r_2| + \frac{\Omega_{ge, p}}{2} \sigma^x_{ge} + \frac{\Omega_{er_1, c}}{2} \sigma^x_{er_1} + \frac{\Omega_{er_2, c}}{2} \sigma^x_{er_2},
$$
where $\sigma^x_{ab}=|a\rangle\langle b| + |b\rangle\langle a|$.


## Valid transitions

At this point it is helpful to identify the non-vanishing transitions in order to limit the size of our atomic Hilbert space. The linearly polarized dipole transitions couple only states for which $\Delta l = \pm 1$, $\Delta j = 0, \pm 1$, and $\Delta m_j = 0$; while the circularly polarized dipole transitions couple states for which $\Delta l = \pm 1$, $\Delta j = 0, \pm 1$, and $\Delta m_j = \pm 1$ (the sign of $\Delta m_j$ corresponding to the handedness of the photon).

Besides those states coupled via the dominant Rabi terms the system may only couple to other states
populated via non-negligible decay processes. Towards an understanding of dominant decay channels
we tabulate below a few of these transition rates. Note that for each given transition listed
there exists also a transition of equal rate with negated initial and final $m_j$ values.

| Upper State, $m_j$ | Lower State, $m_j$ | $\Gamma$ ($\text{us}^{-1}$) |
|-------------|-------------|--------|
| 6P3/2, 3/2  | 6S1/2, 1/2  | 3.3E1  |
| 6P3/2, 1/2  | 6S1/2, 1/2  | 2.2E1  |
| 6 ... 15    | 6 ... 15    | $\mathcal{O}$(3E-1...3E1) |
| $\geq 16$   | 6 ...       | $\mathcal{O}$(<3E-1) |

The behaviour of decay from states with increasing principal number to the ground state is quantified
in the table below,

| Upper State, $m_j$ | Lower State, $m_j$ | $\Gamma$ ($\text{us}^{-1}$) |
|-------------|-------------|--------|
| 6P3/2, 1/2  | 6S1/2, 1/2  | 2.2E1  |
| 16P3/2, 1/2 | 6S1/2, 1/2  | 1.5E-2 |
| 26P3/2, 1/2 | 6S1/2, 1/2  | 2.4E-3 |
| 36P3/2, 1/2 | 6S1/2, 1/2  | 7.9E-4 |

and the behaviour of decay between adjacent energy levels with increasing principal number is
as follows:

| Upper State, $m_j$ | Lower State, $m_j$ | $\Gamma$ ($\text{us}^{-1}$) |
|-------------|-------------|--------|
| 7P3/2, 1/2  | 6S1/2, 1/2  | 1.2E1  |
| 17P3/2, 1/2 | 16S1/2, 1/2 | 9.0E-4 |
| 27P3/2, 1/2 | 26S1/2, 1/2 | 5.1E-5 |
| 37P3/2, 1/2 | 36S1/2, 1/2 | 8.2E-6 |

The magnitudes of these transitions may serve to guide the truncation size of the atomic Hilbert space.


## Spontaneous emission and the master equation

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
\dot{\rho}_{gg} &= i \frac{\Omega_{ge,p}}{2} \left(\rho_{ge} - \rho_{eg}\right) + \Gamma_{ge} \rho_{ee} + \sum_{\alpha>g} \Gamma_{g\alpha} \rho_{\alpha\alpha} \\
\dot{\rho}_{eg;v} &= i \delta_p \rho_{eg;v} + i \frac{\Omega_{ge,p;v}}{2} \left(\rho_{ee} - \rho_{gg}\right) - \frac{\Gamma_{ge}}{2} \rho_{eg} \\
\dot{\rho}_{ee;v} &= i \frac{\Omega_{ge,p;v}}{2} \left(\rho_{eg;v} - \rho_{ge;v}\right) + i \sum_r \frac{\Omega_{er,c;v}}{2}\left(\rho_{er;v}-\rho_{re;v}\right) \nonumber \\
&\quad +\sum_{\alpha>e} \Gamma_{e\alpha} \rho_{\alpha\alpha;v}  - \sum_{\alpha<e} \Gamma_{\alpha e} \rho_{ee;v} \\
\dot{\rho}_{re} &= i \left(\delta_{c,r}-\delta_p\right) \rho_{re} + i \frac{\Omega_{er,c}}{2} \left(\rho_{rr}-\rho_{ee}\right) - \frac{\Gamma_{er}}{2} \rho_{re} \\
\dot{\rho}_{rr} &= i \frac{\Omega_{er,c}}{2}\left(\rho_{re} - \rho_{er} \right) \nonumber \\
&\quad + \sum_{\alpha>r} \Gamma_{r\alpha} \rho_{\alpha\alpha} - \sum_{\alpha<r} \Gamma_{\alpha r} \rho_{rr} \\
\dot{\rho}_{aa} &= \sum_{b>a} \Gamma_{a b} \rho_{bb} - \sum_{b<a} \Gamma_{ba} \rho_{aa}
\end{align}
$$
---
$$
\begin{align}
\dot{\rho}_{rg} &= i\delta_{c,r} \rho_{rg} - \frac{\Gamma_{gr}}{2} \rho_{rg} \\
\dot{\rho}_{a g} &= - \frac{\Gamma_{g a}}{2} \rho_{a g} \\
\dot{\rho}_{a e} &= -i \delta_p \rho_{a e} - \frac{\Gamma_{a e}}{2} \rho_{a e} \\
\dot{\rho}_{a r} &= -i \delta_{c,r} \rho_{a r} - \frac{\Gamma_{a r}}{2} \rho_{a r} \\
\dot{\rho}_{a, b \neq a} &= -\frac{\Gamma_{ab}}{2} \rho_{ab}
\end{align}
$$

where here $r$ denotes either of $r_1,r_2$; $\alpha$ denotes any state; $a,b$ denote any peripheral states; i.e., states outside of $g,e,r_1,r_2$; and an index inequality (e.g., $\alpha>g$) refers to the energy level ordering of the corresponding states. The equations are partitioned by a horizontal dividing line for the reason explained below.


## Complexity of evolution

The above master equations describe the disjoint evolution of two independent collections of density matrix elements. On the one hand, we have the set which includes the primary populations, a subset of the coherences between the primary populations, and the peripheral populations; i.e, the set of $\rho_{\alpha\alpha}, \rho_{ge}, \rho_{e r_1}, \rho_{e r_2}$. These quantities evolve into each other according to the first collection of equations before the horizontal dividing line above.

On the other hand, we have any coherences involving the peripheral states, as well as coherences between ground and Rydberg states; i.e., the set of $\rho_{\alpha a}, \rho_{g r}$. These quantities evolve independently and may be trivially propagated from their initial conditions.

In total, for $N$ peripheral states, we must evolve a set of $N+7$ first-order linear differential equations in $N+7$ unknowns.



## Steady-state solution

Let us denote the system of $N+7$ by $N+7$ equations in as many unkowns $\rho$ by the equation
$$
\dot{\rho} = A \rho.
$$
The steady-state condition is that $0 = A \rho$, which implies one of two cases:
+ $A$ is invertible and thus $\rho$ vanishes, or
+ $A$ is non-invertible and $\rho$ lives in its null space.

In the former case there exists no steady-state solution. In the latter case, if the null space of $A$ is one-dimensional then $\rho$ is determined by the additional constraint that $\text{Tr}[\rho]=1$. If the null space of $A$ is higher than one-dimensional then multiple steady-state solutions exist; however, the initial conditions of the system determine which steady state the system will evolve into.

For our systems we will, in fact, only encounter the case where $A$ has a one-dimensional null space such that $\rho$ is determined by the two equations $A\rho=0$ and $\text{Tr}[\rho]=1$; together, these two equations form $N+8$ equations in $N+7$ unknowns; however, row reduction of $A\rho=0$ will lead to a redundant $0=0$ equation which may be tossed out and replaced with the $\text{Tr}[\rho]=1$ condition. Formally, it is simplest to construct a new $A'$ comprising $A$ augmented with a row defining $\text{Tr}[\rho]$; and $b'$ a vector of zeros augmented with the trace condition 1. Then, we solve the normal equations
$$
{A^\prime}^\text{T} A^\prime \rho = {A^\prime}^\text{T} b^\prime
$$
for the steady density operator state.


## Absorption spectrum

The primary observable of interest is the absorption spectrum of the probe field. Analysis of electromagnetic propagation through a medium finds that the absorption of a field nearly resonant with some transition is equal to the imaginary part of the coherence between the corresponding two states; i.e., for the probe field near the $g,e$ transition, the absorption is given by the imaginary part of $\rho_{eg}$. From the above steady-state equations we may solve for
$$
\rho_{eg} = \frac{\Omega_{ge,p}\left( \rho_{gg} - \rho_{ee} \right)}{2 \delta_p + i \Gamma_{ge}}
$$
such that the absorption $P_\text{abs}$ is
$$
P_\text{abs} = \frac{\Omega_{ge} \Gamma_{ge} / 4}{\delta^2_p + (\Gamma_{ge}/2)^2} \left(\rho_{ee}-\rho_{gg}\right).
$$
The multiplicative factor out front is a Lorentzian in the detuning of the probe field with FWHM equal to $\Gamma_{ge}$, corresponding to a typical line absorption shape. In the case of electromagnetically induced transparency, however, a critical value of $\delta_p$ drives the system into a dark state having no population in $\rho_{ee}$ or in $\rho_{gg}$; the absorption thus vanishes within a very narrow bandwidth where the rightmost multiplier vanishes. THIS NEEDS ANALYSIS...