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
\tilde{\omega}_i =
\begin{cases}
\omega_g, & i = g \\
\omega_g + \omega_p, & i = e \\
\omega_g + \omega_p + \omega_c, & i \in R,
\end{cases}
$$
where $R$ is the manifold of Rydberg states. The Hamiltonian is transformed to
$$
\begin{align}
H^\prime &= U H U^\dagger - i U \frac{dU^\dagger}{dt} \\
&= -\sum_i \delta_i |i\rangle\langle i| + \sum_{k,i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i(\omega_k - \tilde{\Delta}_{ij}) t} |i\rangle\langle j| + \text{h.c.}\right)
\end{align}
$$
where $\tilde{\Delta}_{ij} = \tilde{\omega}_i - \tilde{\omega}_j$.

We now drop fast oscillating terms; i.e., terms for which $\omega_k$ is far from $\tilde{\Delta}_{ij}$. Physically, this means that we retain the interaction between the pump field and the $g,e$ transition; and between the control field and any $e,r\in R$ transitions. Explicitly,
$$
H^\prime = -\delta_p |e\rangle\langle e| + \sum_{r} \left(-\delta_{c,r} |r\rangle\langle r| + \frac{\Omega_{er, p}}{2} \sigma^x_{er}\right).
$$
