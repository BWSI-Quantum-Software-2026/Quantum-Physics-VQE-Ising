# VQE for the Transverse-Field Ising Model — Final Project

## What this project is

We implement the **Variational Quantum Eigensolver (VQE)** algorithm in **Q#** (Microsoft QDK) to find the ground-state energy of the **Transverse-Field Ising Model (TFIM)**, a small quantum spin-chain system used to model magnetic materials. We compare multiple versions of the quantum circuit (different ansatz depths, and open vs. periodic boundary conditions) against each other and against the classically-computed exact answer, to see how circuit design choices trade off accuracy against computational resource cost.

This project is based on and takes inspiration from:
> *Variational Quantum Eigensolver: A Comparative Analysis of Classical and Quantum Optimizer Methods*, [arXiv:2412.19176](https://arxiv.org/abs/2412.19176)

**Course requirement note:** the assignment originally allowed any framework; our course requires the implementation to be done in **Q#** specifically, not Qiskit. This shapes some of the architecture below — Q# doesn't have a built-in Hamiltonian/VQE class the way Qiskit does, so the Hamiltonian is represented as an explicit list of measurable Pauli terms, and the classical optimization loop lives in a Python host program that calls into Q# (see below).

---

## Background: the physics and the algorithm

**The Ising model.** Picture a chain of spins (tiny magnets), each pointing "up" or "down." Neighboring spins influence each other (a coupling), and an external field pushes on each spin individually. The competition between these two effects is what makes the model physically interesting — it's a standard way to study how magnetic order emerges or breaks down.

**The Hamiltonian.** The total energy operator for our system is:

```
H = -J * Σ Zᵢ Zᵢ₊₁  -  h * Σ Xᵢ
```

- `J` — nearest-neighbor coupling strength (how strongly neighboring spins want to align)
- `h` — transverse-field strength (how strongly each spin is pushed by the external field)
- `Z`, `X` — Pauli operators, representing measurement/interaction along two different (perpendicular) axes of a spin

Finding the **ground state** (lowest-energy configuration) of this Hamiltonian is the target problem. For small systems it can be solved exactly by classical diagonalization; the point of VQE is to show a quantum-circuit-based method can approximate the same answer, in a way that in principle scales better than classical diagonalization as the system grows.

**VQE, briefly.** VQE is a hybrid quantum-classical algorithm:
1. A parameterized quantum circuit (the "ansatz") prepares a trial state.
2. The quantum computer measures the expectation value of the Hamiltonian on that state (the "energy" of the trial state).
3. A classical optimizer adjusts the circuit's parameters to try to lower that energy.
4. Repeat until the energy stops improving — the final value approximates the ground-state energy.

---

## Project structure

```
Ansatz.qs               Parameterized quantum circuit (the VQE trial-state preparation)
Hamiltonian.qs           Hamiltonian construction + energy-measurement machinery
VQE.qs                   Wires the ansatz and Hamiltonian together into one energy function
estimate_resources.py    Python driver: calls the Q# Resource Estimator on the VQE circuit
exact.py                 Classical baseline: exact diagonalization of the Hamiltonian matrix
requirements.txt         Python dependencies for the host-side scripts
run_vqe.py               Python driver: runs the full VQE optimization loop (scipy + Q#)
test_exact.py            Unit tests (pytest) for exact.py
```

Q# and Python are split by role: **Q# handles everything quantum** (state prep, measurement), and **Python handles everything classical** (numerical optimization, exact diagonalization, resource-estimation reporting), calling into Q# through the `qdk.qsharp` package.

---


### `Hamiltonian.qs` — the Hamiltonian and energy measurement

This file defines the physical system that the Variational Quantum Eigensolver (VQE) is trying to solve. In VQE, the goal is to find the lowest-energy (ground) state of a Hamiltonian, so this file provides the mathematical description of the Transverse-Field Ising Model (TFIM). It is also responsible for measuring the energy of a quantum state prepared by the ansatz. Since VQE repeatedly evaluates the Hamiltonian throughout the optimization process, this file serves as the foundation of the entire quantum simulation.

- **`HamiltonianTerm`** — a struct bundling one weighted Pauli term of the Hamiltonian (coefficient, which Pauli operators, which qubits they act on).
- **`GenerateTFIMHamiltonian(numQubits, J, h, periodic)`** — builds the full list of Hamiltonian terms for a spin chain of any size.
  - `periodic = false` → open chain (free ends), `numQubits - 1` coupling bonds.
  - `periodic = true` → periodic ring (matches the reference paper's topology), adds one extra "wrap-around" bond connecting the last spin back to the first. (Guarded so a 2-qubit chain doesn't double-count its only bond as both the open bond and the wrap bond.)
- **`SelectTargets`** — helper that picks out the specific qubits a term acts on.
- **`MeasureTermEigenvalue`** — measures one Hamiltonian term via Q#'s built-in joint Pauli measurement, returning ±1.
- **`EstimateEnergy(prepareState, numQubits, terms, shotsPerTerm)`** — the core energy-evaluation function: given *any* state-preparation operation, repeatedly prepares the state and measures every Hamiltonian term to estimate `⟨H⟩`.
- **`PrepareAllZero`** — a trivial reference state used only for verification (its energy is analytically known by hand).
- **`Main()`** — a demo entry point that builds both the open-chain and periodic-ring Hamiltonians and prints their energies side by side on the reference state.
- **7 unit tests** (`@Test("QuantumSimulator")` + `Fact`): term-count correctness for both boundary conditions, the 2-qubit edge case, wrap-bond correctness, field-term correctness, and — most importantly — energy estimates matching the analytically-known exact value for both boundary conditions (this validates the physics, not just the code structure).

### `Ansatz.qs` — the parameterized trial-state circuit

This file implements the parameterized quantum circuit, or **ansatz**, used by the Variational Quantum Eigensolver. The ansatz prepares a trial quantum state whose parameters are adjusted by the classical optimizer until the measured energy is minimized. Because the quality of the ansatz directly affects the accuracy of VQE, this file is responsible for generating the quantum states that are evaluated throughout the optimization process.


- **`RealAmplitudesAnsatz(theta, qubits, reps)`** — a hardware-efficient ansatz: `reps` repeated layers, each layer applying one tunable `Ry` rotation per qubit, followed by a chain of `CNOT`s linking each qubit to the next. `theta` holds `reps * numQubits` angles total (one per qubit per layer). Marked `is Adj + Ctl` so it composes cleanly with other Q# operations that need adjoint/controlled variants.
- This is the "how deep is the circuit" knob for the project's comparison: varying `reps` directly varies circuit depth, gate count, and (expected) accuracy — the sweep the resource-estimation and results sections are built around.


### `VQE.qs` — connecting ansatz to Hamiltonian

This file combines the ansatz and Hamiltonian into a single quantum routine that can be called by the classical optimizer. Rather than performing the optimization itself, it prepares a parameterized quantum state, measures its energy with respect to the Hamiltonian, and returns that energy as a single value. During VQE, this operation is called hundreds or even thousands of times as the optimizer searches for the lowest-energy state.


- **`VQEEnergy(theta, numQubits, J, h, reps, periodic, shotsPerTerm)`** — builds the Hamiltonian, builds the ansatz-prepared state from `theta`, and returns the estimated energy via `EstimateEnergy`. This is the single function the classical optimizer calls, over and over, with different `theta` each time.

### `exact.py` — classical ground-truth baseline

VQE is an approximate algorithm, so its results need to be compared against an exact solution whenever possible. This file provides that classical benchmark by constructing the same Hamiltonian as a NumPy matrix and solving it using exact diagonalization. The resulting ground-state energy is used to evaluate how accurately the quantum algorithm performs.


- **`to_matrix(H)`** — normalizes a Hamiltonian (however it's represented) to a plain NumPy array.
- **`exact_ground_energy(H)`** — lowest eigenvalue via `np.linalg.eigvalsh`.
- **`exact_ground_state(H)`** — lowest eigenvalue *and* its eigenvector.
- **`exact_spectrum(H)`** — the full sorted list of energy levels, not just the ground state.
- Includes a runnable self-check (`__main__` block) verifying two known cases by hand: `H = -X` (ground energy `-1`) and the 2-spin Ising Hamiltonian `H = -Z₀Z₁ - X₀ - X₁` (ground energy `-√5 ≈ -2.236`), which is the textbook sanity-check value for this model.

### `test_exact.py` — unit tests for the classical baseline

This file contains the Python unit tests for the classical solver. Since the exact solution is used as the benchmark throughout the project, it is important to verify that the classical implementation is producing the correct results before comparing it against the quantum algorithm.


6 pytest tests: single-spin `-X` and `Z` ground energies, the known `-√5` two-spin value, spectrum length/ordering (2 qubits → 4 levels, sorted ascending), the eigen-equation check (`H·ψ = E·ψ` for the returned ground state), and a trivial diagonal-matrix case as a basic correctness anchor.

### `run_vqe.py` — the full VQE optimization loop (Python host + Q#)

This file is the main driver for the project and implements the hybrid quantum-classical optimization loop. It connects the Python optimizer to the Q# quantum routines by repeatedly proposing new ansatz parameters, calling the quantum energy evaluation, and updating the parameters until the estimated energy converges. This script is responsible for executing the complete VQE algorithm from start to finish.


- **`load_qsharp()`** — compiles `Ansatz.qs`, `Hamiltonian.qs`, `VQE.qs` into the running Python session via `qdk.qsharp`.
- **`energy(theta, n, J, h, reps, periodic, shots)`** — thin wrapper that formats `theta` into a Q# call to `VQEEnergy` and returns the result.
- **`one_run(...)`** — a single optimization attempt: `scipy.optimize.minimize` with `COBYLA` proposes angles, `energy()` scores them, and the best angles found are re-measured with more shots at the end (since the noisy sampled minimum during optimization is statistically biased low, so it needs a higher-precision remeasurement to report honestly).
- **`run_vqe(...)`** — wraps `one_run` with **multiple random restarts** (default 5), since COBYLA can get stuck in a local minimum depending on its starting angles; picks the best restart, then does one final high-shot remeasurement of the winning angles for the reported number.
- **`at(op, i, n)`** and **`tfim_matrix(n, J, h, periodic)`** — build the *same* Hamiltonian as `Hamiltonian.qs`, but as a plain NumPy matrix, so `exact_ground_energy` can be run on it for comparison. (Worth double-checking these two Hamiltonian constructions — the Q# one and this NumPy one — stay in sync if either changes.)
- **`__main__` block** — runs VQE for `n = 2, 3`, compares against the exact floor, and prints a table of VQE energy, exact energy, absolute error, and number of energy evaluations used.

### `estimate_resources.py` — resource estimation

This file analyzes the quantum resources required to run the VQE circuit using the Microsoft Quantum Development Kit Resource Estimator. Rather than evaluating the accuracy of the algorithm, it estimates how expensive the quantum circuit would be on a fault-tolerant quantum computer by reporting gate counts, qubit counts, runtime estimates, and other hardware metrics. These results are used to study how the cost of the algorithm scales as the problem size increases.


- **`load_qsharp(filename)`** — compiles a single `.qs` file so its operations are callable.
- **`estimate_resources(entry_expr, params=None)`** — calls Q#'s built-in Resource Estimator on a Q# call expression (e.g. `"VQEEnergy(...)"`) and pulls out: logical qubit count, rotation-gate count, T-gate count, measurement count, physical qubit count, runtime, code distance, and number of T-factories.
- **`print_estimate(...)`** — pretty-prints one estimate.
- **`__main__` block** — sweeps system size `n ∈ {2, 3, 4}` and ansatz depth `reps ∈ {1, 2, 3}`, running the resource estimator on `VQEEnergy` at each combination. This is exactly the sweep the project's "circuit depth vs. resource cost" comparison needs.

### `requirements.txt`

This file lists all of the Python dependencies required to run the classical portion of the project. Installing these packages ensures that the optimization routines, testing framework, numerical computations, resource estimation scripts, and the Python–Q# interface all function correctly without additional manual setup.

Lists the Python dependencies needed to run the host-side scripts (the `qdk.qsharp` package, `numpy`, `scipy`, `pytest`, based on what `run_vqe.py`, `exact.py`, `estimate_resources.py`, and `test_exact.py` import).


