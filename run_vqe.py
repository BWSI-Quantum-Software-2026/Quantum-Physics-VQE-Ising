# python driver, scipy > runs Q# measures the energy

from pathlib import Path

import numpy as np
import qdk.qsharp as qs
from scipy.optimize import minimize
import matplotlib.pyplot as plt

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


def one_run(n, J, h, shots, maxiter, start, check_shots):
    # scipy picks angles then we have Q# to score them, finally repeat until the energy stops dropping
    history = []

    def objective(theta):
        value = energy(theta, n, J, h, shots)
        history.append(value)
        print(len(history), value)
        return value

    result = minimize(objective, start, method="COBYLA",
                      options={"maxiter": maxiter})

    # the best sample seen is biased low by noise, so remeasure the final angles
    return energy(result.x, n, J, h, check_shots), result.x, history


def run_vqe(n, J=1.0, h=1.0, shots=4000, maxiter=300, seed=0,
            restarts=5, check_shots=20000, final_shots=100000):
    # COBYLA can settle in a local minimum, so we try several random starts
    rng = np.random.default_rng(seed)
    runs = [one_run(n, J, h, shots, maxiter,
                    rng.uniform(0, 2 * np.pi, 2 * n), check_shots)
            for _ in range(restarts)]

    _, theta, history = min(runs, key=lambda r: r[0])

    # measure the winner once more so the reported number is not the pick of a batch
    return energy(theta, n, J, h, final_shots), theta, history


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

    # store results for graphs
    n_values = []
    vqe_values = []
    exact_values = []
    errors = []
    eval_counts = []

    for n in (2, 3, 4):
        best, theta, history = run_vqe(n)
        floor = exact_ground_energy(tfim_matrix(n, 1.0, 1.0))
        print(f"{n:>2} {best:>9.4f} {floor:>9.4f} {abs(best - floor):>8.4f} {len(history):>6}")

        #save values for graphs
        n_values.append(n)
        vqe_values.append(best)
        exact_values.append(floor)
        errors.append(abs(best-floor))
        eval_counts.append(len(history))

        # ---------- Convergence Plot ----------
        plt.figure(figsize=(7,5))

        plt.plot(
            range(1, len(history)+1),
            history,
            marker='o',
            linewidth=2,
            label="VQE Energy"
        )

        plt.axhline(
            y=floor,
            color='red',
            linestyle='--',
            linewidth=2,
            label="Exact Energy"
        )

        plt.text(
            len(history),
            history[-1],
            f"Final = {history[-1]:.4f}",
            fontsize=9,
            ha="right",
            va="bottom"
        )

        plt.text(
            len(history),
            floor,
            f"Exact = {floor:.4f}",
            fontsize=9,
            ha="right",
            va="top",
            color="red"
        )

        plt.xlabel("Optimization Iteration")
        plt.ylabel("Energy")
        plt.title(f"VQE Convergence (n={n})")

        plt.legend()
        plt.grid(True)

        plt.tight_layout()

        plt.savefig(f"convergence_n{n}.png", dpi=300)

        plt.show()
