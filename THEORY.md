# THEORY -- MAXWELL-BLOCH EQUATIONS FOR A MULTI-LEVEL RYDBERG ATOM

## Rotating frame and associated approximations

We generalise here the rotating frame method to a system of $N$ atomic levels having energies $\omega_i$ as dressed by fields described in README. Let the Rabi frequency of the $k\text{th}$ field coupling $|i\rangle$ and $|j\rangle$ be $\Omega_{ij,k}$. Letting $\hbar\rightarrow 1$ the total Hamiltonian is
$$
H = \sum_{i=1}^N \omega_i |i\rangle \langle i| + \sum_k \sum_{i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i \omega_k t} |i\rangle\langle j| + \text{h.c.}\right).
$$
Note that a subscript $i$ refers to an atomic energy level, while a subscript $k$ refers to a laser field angular frequency.

We construct a rotating frame unitary transformation which strategically factors off driving laser angular frequencies from the various atomic transitions. Let
$$
|\psi(t)\rangle \rightarrow |\psi^\prime(t)\rangle = U(t) |\psi(t)\rangle,
$$
where
$$
U(t) = \sum_{i=1}^N e^{i \tilde{\omega}_i t} | i \rangle\langle i |
$$
for $\tilde{\omega}_i$ a set of rotating frame angular frequencies defined by
$$
\begin{align}
\tilde{\omega}_g &= \omega_g \\
\tilde{\omega}_e &= \omega_g + \omega_p \\
\tilde{\omega}_r &= \omega_e + \omega_c, r \in R
\end{align}
$$
where $R$ is the manifold of Rydberg states. The Hamiltonian is transformed to
$$
\begin{align}
H^\prime &= U H U^\dagger - i U \frac{dU^\dagger}{dt} \\
&= -\sum_i \delta_i |i\rangle\langle i| + \sum_{k,i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i(\omega_k - \tilde{\Delta}_{ij}) t} |i\rangle\langle j| + \text{h.c.}\right)
\end{align}
$$
where $\tilde{\Delta}_{ij} = \tilde{\omega}_i - \tilde{\omega}_j$.

We now drop fast oscillating terms; i.e., terms for which $\omega_k$ is far from $\tilde{\Delta}_{ij}$. Physically, this means that we retain the interaction between the pump field and the $g,e$ transition; and between the control field and any $e,r\in R$ transitions, such that
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
\dot{\rho}_{gg} &= i \frac{\Omega_{ge,p}}{2} \left(\rho_{ge} - \rho_{eg}\right) + \Gamma_{eg} \rho_{ee} + \sum_r \Gamma_{rg} \rho_{rr} \\
\dot{\rho}_{eg} &= i \delta_p \rho_{eg} + i \frac{\Omega_{ge,p}}{2} \left(\rho_{ee} - \rho_{gg}\right) - \frac{\Gamma_{eg}}{2} \rho_{eg} \\
\dot{\rho}_{rg} &= \\
\dot{\rho}_{ee} &= i \frac{\Omega_{ge,p}}{2} \left(\rho_{eg} - \rho_{ge}\right) + i \sum_r \frac{\Omega_{er,c}}{2}\left(\rho_{er}-\rho_{re}\right) - \Gamma_{eg} \rho_{ee} + \sum_r \Gamma_{re} \rho_{rr}
\end{align} \\

$$