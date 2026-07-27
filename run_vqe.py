# python driver, scipy > runs Q# measures the energy

from pathlib import Path

import numpy as np
import qdk.qsharp as qs
from scipy.optimize import minimize

from exact import exact_ground_energy

QS_FILES = ("Ansatz.qs", "Hamiltonian.qs", "VQE.qs")

I = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def load_qsharp():
    # compile the Q# files so VQEEnergy can be called
    folder = Path(__file__).parent
    for name in QS_FILES:
        qs.eval((folder / name).read_text(encoding="utf-8"))


def energy(theta, n, J, h, shots):
    # ask Q# for the energy of the state these angles build
    angles = ", ".join(str(float(t)) for t in theta)
    return qs.eval(f"VQEEnergy([{angles}], {n}, {J}, {h}, {shots})")


def run_vqe(n, J=1.0, h=1.0, shots=8000, maxiter=300, seed=0):
    # scipy picks angles, Q# scores them, repeat until the energy stops dropping
    rng = np.random.default_rng(seed)
    start = rng.uniform(0, 2 * np.pi, 2 * n)
    history = []

    def objective(theta):
        value = energy(theta, n, J, h, shots)
        history.append(value)
        return value

    result = minimize(objective, start, method="COBYLA",
                      options={"maxiter": maxiter})
    return result.fun, result.x, history


def at(op, i, n):
    # place a single qubit operator on spin i of an n spin chain
    m = np.array([[1]], dtype=complex)
    for k in range(n):
        m = np.kron(m, op if k == i else I)
    return m


def tfim_matrix(n, J, h):
    # same Hamiltonian as Hamiltonian.qs but as a numpy matrix
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(n - 1):
        H += -J * (at(Z, i, n) @ at(Z, i + 1, n))
    for i in range(n):
        H += -h * at(X, i, n)
    return H


if __name__ == "__main__":
    load_qsharp()

    print(f"{'n':>2} {'VQE':>9} {'exact':>9} {'error':>8} {'evals':>6}")
    for n in (2, 3):
        best, theta, history = run_vqe(n)
        floor = exact_ground_energy(tfim_matrix(n, 1.0, 1.0))
        print(f"{n:>2} {best:>9.4f} {floor:>9.4f} {abs(best - floor):>8.4f} {len(history):>6}")
