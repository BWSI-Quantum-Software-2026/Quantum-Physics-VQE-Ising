# sweeps the magnetic field to look for the quantum phase transition

import json
import time
from pathlib import Path

import numpy as np
import qdk.qsharp as qs

from exact import exact_ground_energy
from run_vqe import energy, nft, tfim_matrix

QS_FILES = ("Ansatz.qs", "Hamiltonian.qs", "VQE.qs", "Observables.qs")

N = 4
REPS = 3
PERIODIC = True
J = 1.0

# sweep downward, the strong field end is easy so each point warm starts the next
FIELDS = np.round(np.arange(2.0, -0.001, -0.1), 3)

SHOTS = 3000
MAXITER = 200
OBS_SHOTS = 20000
OUT = Path(__file__).parent / "sweep.json"


def load_qsharp():
    folder = Path(__file__).parent
    for name in QS_FILES:
        qs.eval((folder / name).read_text(encoding="utf-8"))


def observables(theta, n, reps, periodic, shots):
    # average transverse magnetisation and neighbour correlation
    angles = ", ".join(str(float(t)) for t in theta)
    flag = "true" if periodic else "false"
    mx = qs.eval(f"MeasureTransverse([{angles}], {n}, {reps}, {shots})")
    zz = qs.eval(f"MeasureCoupling([{angles}], {n}, {reps}, {flag}, {shots})")
    return mx, zz


def sweep():
    load_qsharp()
    rng = np.random.default_rng(0)
    theta = rng.uniform(0, 2 * np.pi, REPS * N)
    records = []

    for h in FIELDS:
        started = time.time()

        def objective(t):
            return energy(t, N, J, h, REPS, PERIODIC, SHOTS)

        theta = nft(objective, theta, MAXITER)

        vqe = energy(theta, N, J, h, REPS, PERIODIC, OBS_SHOTS)
        mx, zz = observables(theta, N, REPS, PERIODIC, OBS_SHOTS)
        floor = exact_ground_energy(tfim_matrix(N, J, h, PERIODIC))

        records.append({
            "n": N,
            "reps": REPS,
            "periodic": PERIODIC,
            "J": J,
            "h": float(h),
            "vqe_energy": vqe,
            "exact_energy": floor,
            "error": abs(vqe - floor),
            "transverse_magnetisation": mx,
            "neighbour_correlation": zz,
            "angles": list(theta),
            "seconds": time.time() - started,
        })

        OUT.write_text(json.dumps(records, indent=2), encoding="utf-8")
        print(f"h={h:>5.2f}  E/n {vqe / N:>8.4f}  exact {floor / N:>8.4f}  "
              f"<X> {mx:>6.3f}  <ZZ> {zz:>6.3f}  ({time.time() - started:.0f}s)")

    return records


if __name__ == "__main__":
    runs = sweep()
    print(f"\nsaved {len(runs)} field points to {OUT.name}")
