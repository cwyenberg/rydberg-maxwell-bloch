# THEORY -- MAXWELL-BLOCH EQUATIONS FOR A MULTI-LEVEL RYDBERG ATOM

## Rotating frame and associated approximations

We generalize here the rotating frame method to a system of $N$ atomic levels having energies $\omega_i$ and dressed by one ore more laser fields. Let the Rabi frequency of the $k$th field coupling $|i\rangle$ and $|j\rangle$ be $\Omega_{ij,k}$ with detuning $\Delta_k = \omega_{ij} - \omega_k$ for $\omega_{ij} = \omega_i - \omega_j$. Letting $\hbar\rightarrow 1$, the total Hamiltonian is
$$
H = \sum_{i=1}^N \omega_i |i\rangle \langle i| + \sum_k \sum_{i<j} \left(\frac{\Omega_{ij,k}}{2} e^{-i \omega_k t} |i\rangle\langle j| + \text{h.c.}\right).
$$
Note that a subscript $i$ refers to an atomic energy level, while a subscript $k$ refers to a laser field angular frequency.

The rotating frame unitary transformation is defined as
$$
U(t) = \sum_{i=1}^N e^{-i \omega_i t} | i \rangle\langle i |,
$$
such that
$$
|\psi(t)\rangle \rightarrow |\psi^\prime(t)\rangle = U(t) |\psi(t)\rangle,
$$
and the Hamiltonian transforms to
$$
H^\prime = U H U^\dagger - i U \frac{dU^\dagger}{dt}.
$$
