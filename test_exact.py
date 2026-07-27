# unit tests for exact.py

import numpy as np
import pytest

from exact import (
    exact_ground_energy,
    exact_ground_state,
    exact_spectrum,
)

# Single-qubit Pauli matrices, used to build test Hamiltonians.
I = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def two_spin_ising():
    # H = -Z0Z1 - X0 - X1 (J = h = 1). Known ground energy is -sqrt(5).
    return -np.kron(Z, Z) - np.kron(I, X) - np.kron(X, I)


def test_single_spin_x_field():
    # H = -X has levels -1 and +1, so ground energy is -1.
    assert exact_ground_energy(-X) == pytest.approx(-1.0)


def test_single_spin_z():
    # H = Z has levels -1 and +1, so ground energy is -1.
    assert exact_ground_energy(Z) == pytest.approx(-1.0)


def test_two_spin_ising_known_value():
    assert exact_ground_energy(two_spin_ising()) == pytest.approx(-np.sqrt(5))


def test_spectrum_is_sorted_and_right_length():
    # 2 qubits means 2^2 = 4 energy levels, sorted low to high.
    levels = exact_spectrum(two_spin_ising())
    assert len(levels) == 4
    assert np.all(np.diff(levels) >= 0)


def test_ground_state_solves_eigen_equation():
    # The returned state psi must satisfy H @ psi = E * psi.
    H = two_spin_ising()
    energy, psi = exact_ground_state(H)
    assert np.allclose(H @ psi, energy * psi)
    assert energy == pytest.approx(exact_ground_energy(H))


def test_diagonal_matrix_known_eigenvalues():
    # A diagonal matrix's eigenvalues are just its diagonal entries.
    H = np.diag([3.0, -2.0, 5.0, 0.0])
    assert exact_ground_energy(H) == pytest.approx(-2.0)
