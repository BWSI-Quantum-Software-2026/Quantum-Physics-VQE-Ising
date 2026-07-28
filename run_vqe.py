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


def energy(theta, n, J, h, reps, periodic, shots):
    # ask Q# for the energy of the state these angles build
    angles = ", ".join(str(float(t)) for t in theta)
    flag = "true" if periodic else "false"
    return qs.eval(f"VQEEnergy([{angles}], {n}, {J}, {h}, {reps}, {flag}, {shots})")


def nft(objective, start, maxiter):
    # with only Ry gates the energy is an exact sine wave in each angle,
    # so three measurements pin that wave down and we jump to its minimum
    theta = [float(t) for t in start]
    sweeps = max(1, maxiter // (3 * len(theta)))

    for _ in range(sweeps):
        for d in range(len(theta)):
            theta[d] = 0.0
            m0 = objective(theta)
            theta[d] = np.pi / 2
            mp = objective(theta)
            theta[d] = -np.pi / 2
            mm = objective(theta)
            theta[d] = -np.pi / 2 - np.arctan2(2 * m0 - mp - mm, mp - mm)

    return np.array(theta)


def one_run(n, J, h, reps, periodic, shots, maxiter, start, check_shots,
            method="COBYLA"):
    # the optimizer picks angles then we have Q# to score them, repeat until the energy stops dropping
    history = []

    def objective(theta):
        value = energy(theta, n, J, h, reps, periodic, shots)
        history.append(value)
        return value

    if method == "NFT":
        best = nft(objective, start, maxiter)
    else:
        best = minimize(objective, start, method=method,
                        options={"maxiter": maxiter}).x

    # the best sample seen is biased low by noise, so remeasure the final angles
    return energy(best, n, J, h, reps, periodic, check_shots), best, history


def run_vqe(n, J=1.0, h=1.0, reps=2, periodic=False, shots=4000, maxiter=300,
            seed=0, restarts=5, check_shots=20000, final_shots=100000,
            method="COBYLA"):
    # the optimizer can settle in a local minimum, so we try several random starts
    rng = np.random.default_rng(seed)
    runs = [one_run(n, J, h, reps, periodic, shots, maxiter,
                    rng.uniform(0, 2 * np.pi, reps * n), check_shots, method)
            for _ in range(restarts)]

    _, theta, history = min(runs, key=lambda r: r[0])

    # measure the winner once more so the reported number is not the pick of a batch
    return energy(theta, n, J, h, reps, periodic, final_shots), theta, history


def at(op, i, n):
    # place a single qubit operator on spin i of an n spin chain
    m = np.array([[1]], dtype=complex)
    for k in range(n):
        m = np.kron(m, op if k == i else I)
    return m


def tfim_matrix(n, J, h, periodic=False):
    # same Hamiltonian as Hamiltonian.qs but as a numpy matrix
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(n - 1):
        H += -J * (at(Z, i, n) @ at(Z, i + 1, n))
    if periodic and n > 2:
        H += -J * (at(Z, n - 1, n) @ at(Z, 0, n))
    for i in range(n):
        H += -h * at(X, i, n)
    return H


if __name__ == "__main__":
    load_qsharp()

    print(f"{'n':>2} {'VQE':>9} {'exact':>9} {'error':>8} {'per site':>9} {'evals':>6}")

    # store results for graphs
    n_values = []
    vqe_values = []
    exact_values = []
    errors = []
    eval_counts = []

    for n in (2, 3, 4):
        best, theta, history = run_vqe(n, method="NFT")
        floor = exact_ground_energy(tfim_matrix(n, 1.0, 1.0))
        # energy grows with the chain, so the per site error is the fair one
        print(f"{n:>2} {best:>9.4f} {floor:>9.4f} {abs(best - floor):>8.4f} "
              f"{abs(best - floor) / n:>9.4f} {len(history):>6}")

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

        # keep this one on the left so it does not sit on top of the label above
        plt.text(
            1,
            floor,
            f"Exact = {floor:.4f}",
            fontsize=9,
            ha="left",
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

        # no plt.show, the window blocks the run and the file is already saved
        plt.close()
