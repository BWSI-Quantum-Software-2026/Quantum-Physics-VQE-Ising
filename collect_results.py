# runs the sweep once and saves every number the plots and tables need

import json
import time
from pathlib import Path

from estimate_resources import estimate_resources
from exact import exact_ground_energy
from run_vqe import load_qsharp, run_vqe, tfim_matrix

SIZES = (2, 3, 4)
DEPTHS = (1, 2, 3)
J, H = 1.0, 1.0
PERIODIC = False

SHOTS = 4000
MAXITER = 250
RESTARTS = 4

# run both so the saved data can compare them
METHODS = ("NFT", "COBYLA")

OUT = Path(__file__).parent / "results.json"


def collect():
    load_qsharp()
    records = []

    for n in SIZES:
        floor = exact_ground_energy(tfim_matrix(n, J, H, PERIODIC))

        for reps in DEPTHS:
            theta = ", ".join(["0.5"] * (reps * n))
            flag = "true" if PERIODIC else "false"
            cost = estimate_resources(
                f"VQEEnergy([{theta}], {n}, {J}, {H}, {reps}, {flag}, 1)")

            for method in METHODS:
                started = time.time()
                vqe, angles, history = run_vqe(
                    n, J=J, h=H, reps=reps, periodic=PERIODIC,
                    shots=SHOTS, maxiter=MAXITER, restarts=RESTARTS,
                    method=method)
                elapsed = time.time() - started

                records.append({
                    "n": n,
                    "reps": reps,
                    "method": method,
                    "periodic": PERIODIC,
                    "J": J,
                    "h": H,
                    "vqe_energy": vqe,
                    "exact_energy": floor,
                    "error": abs(vqe - floor),
                    "error_per_site": abs(vqe - floor) / n,
                    "num_parameters": reps * n,
                    "angles": list(angles),
                    "history": history,
                    "seconds": elapsed,
                    "logical_qubits": cost["logical_qubits"],
                    "rotation_gates": cost["rotation_count"],
                    "t_gates": cost["t_count"],
                    "measurements": cost["measurement_count"],
                    "physical_qubits": cost["physical_qubits"],
                    "runtime": cost["runtime_pretty"],
                })

                # save after every run so a crash does not lose the earlier ones
                OUT.write_text(json.dumps(records, indent=2), encoding="utf-8")
                print(f"n={n} reps={reps} {method:>6}  vqe {vqe:.4f}  "
                      f"exact {floor:.4f}  error {abs(vqe - floor):.4f}  ({elapsed:.0f}s)")

    return records


if __name__ == "__main__":
    runs = collect()
    print(f"\nsaved {len(runs)} runs to {OUT.name}")
