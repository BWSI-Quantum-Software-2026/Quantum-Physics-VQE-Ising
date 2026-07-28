# classical baseline, the  ground-state energy from a Hamiltonian matrix

import numpy as np


def to_matrix(H):
    """Return H as a numpy matrix."""
    if hasattr(H, "to_matrix"):
        return np.asarray(H.to_matrix())
    return np.asarray(H)


def exact_ground_energy(H):
    """Lowest energy level """
    # eigvalsh returns real eigenvalues sorted low -> high. [0] is the lowest.
    return float(np.linalg.eigvalsh(to_matrix(H))[0])


def exact_ground_state(H):
    """Return (lowest energy, its state vector)."""
    energies, states = np.linalg.eigh(to_matrix(H))
    return float(energies[0]), states[:, 0]


def exact_spectrum(H):
    """All energy levels, sorted low -> high."""
    return np.linalg.eigvalsh(to_matrix(H))


def tfim_open_chain_energy(n, J, h):
    """Ground energy of the open chain without building the 2^n matrix."""
    # free fermion result, the energy is minus the singular values of this
    # small n by n matrix, so this stays cheap for hundreds of spins
    M = np.diag([float(h)] * n) + np.diag([float(J)] * (n - 1), 1)
    return -float(np.linalg.svd(M, compute_uv=False).sum())


if __name__ == "__main__":
    # Build example Hamiltonians as plain numpy matrices
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # one spin, H = -X, known answer levels are -1 and +1, ground = -1
    H1 = -X
    print("H = -X   levels:", exact_spectrum(H1), " ground:", exact_ground_energy(H1))

    # two spins, H = -Z0Z1 - X0 - X1, known answer  = -sqrt(5)
    H2 = -np.kron(Z, Z) - np.kron(I, X) - np.kron(X, I)
    print("H = -Z0Z1 - X0 - X1   ground:", round(exact_ground_energy(H2), 6),
          " (expected -2.236068)")
