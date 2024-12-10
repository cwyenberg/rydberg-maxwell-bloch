# THEORY -- MAXWELL-BLOCH EQUATIONS FOR A MULTI-LEVEL RYDBERG ATOM

## Rotating frame and associated approximations

We generalize here the rotating frame method to a system of $N$ atomic levels having energies $\omega_i$ and dressed by one ore more laser fields. Let the Rabi frequency of the $k$th field coupling $|i\rangle$ and $|j\rangle$ be $\Omega_{ij,k}$ with detuning $\Delta_{ij,k} = \omega_{ij} - \omega_k$ for $\omega_{ij} = \omega_i - \omega_j$. Letting $\hbar\rightarrow 1$, the total Hamiltonian is
$$
H = \sum_{i=1}^N \omega_i |i\rangle \langle i| + \sum_k \sum_{i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i \omega_k t} |i\rangle\langle j| + \text{h.c.}\right).
$$
Note that a subscript $i$ refers to an atomic energy level, while a subscript $k$ refers to a laser field angular frequency. Herein we will suppress the ranges on these indices.

Let us define a rotating frame unitary transformation which factors off the natural non-interacting rotations of the various atomic energy levels. To that end, let
$$
|\psi(t)\rangle \rightarrow |\psi^\prime(t)\rangle = U(t) |\psi(t)\rangle,
$$
where the rotating frame unitary transformation is defined as
$$
U(t) = \sum_{i=1}^N e^{i \omega_i t} | i \rangle\langle i |,
$$
such that the Hamiltonian transforms to
$$
\begin{align}
H^\prime &= U H U^\dagger - i U \frac{dU^\dagger}{dt} \\
&= \sum_{k,i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i(\omega_k - \Delta_{ik}) t} |i\rangle\langle j| + \text{h.c.}\right)
\end{align}
$$
where $\Delta_{ij} = \omega_i - \omega_j$. Note that the internal atomic energies have vanished in the rotating frame.
