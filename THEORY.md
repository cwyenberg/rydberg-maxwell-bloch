# THEORY -- MAXWELL-BLOCH EQUATIONS FOR A MULTI-LEVEL RYDBERG ATOM

## Rotating frame and associated approximations

In this document, please refer to the README.md file for any undefined quantities.

We generalise here the rotating frame method to a system of $N$ atomic levels having energies $\omega_i$ and dressed by one ore more laser fields. Let the Rabi frequency of the $k\text{th}$ field coupling $|i\rangle$ and $|j\rangle$ be $\Omega_{ij,k}$. Letting $\hbar\rightarrow 1$ the total Hamiltonian is
$$
H = \sum_{i=1}^N \omega_i |i\rangle \langle i| + \sum_k \sum_{i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i \omega_k t} |i\rangle\langle j| + \text{h.c.}\right).
$$
Note that a subscript $i$ refers to an atomic energy level, while a subscript $k$ refers to a laser field angular frequency.

It will be important in what follows to differentiate between three classes of atomic levels: the first, those near the ground energy level, which we denote as belonging to the set $G$; the second, those near the excited level belonging to the set $E$; and the third, those near the rydberg level belonging to the set $R$. This partitioning will define the rotating frame factors associated to each state below.

Let us define a rotating frame unitary transformation which strategically factors off driving laser angular frequencies from the various atomic energy levels. Let
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
\omega_g, & i \in G \\
\omega_g + \omega_p, & i \in E \\
\omega_e + \omega_c, & i \in R
\end{cases}.
$$



such that the Hamiltonian transforms to

$$
\begin{align}
H^\prime &= U H U^\dagger - i U \frac{dU^\dagger}{dt} \\
&= \sum_{k,i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i(\omega_k - \Delta_{ik}) t} |i\rangle\langle j| + \text{h.c.}\right)
\end{align}
$$
where $\Delta_{ij} = \omega_i - \omega_j$. Note that the internal atomic energies have vanished in the rotating frame.
