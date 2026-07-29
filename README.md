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
Ansatz.qs                Parameterized quantum circuit (the VQE trial-state preparation)
Hamiltonian.qs           Hamiltonian construction + energy-measurement machinery
Observables.qs           Measures the physics quantities we plot (magnetisation, correlation)
VQE.qs                   Wires the ansatz and Hamiltonian together into one energy function
collect_results.py       Python driver: runs the size/depth sweep once and saves results.json
estimate_resources.py    Python driver: calls the Q# Resource Estimator on the VQE circuit
exact.py                 Classical baseline: exact diagonalization of the Hamiltonian matrix
requirements.txt         Python dependencies for the host-side scripts
run_vqe.py               Python driver: runs the full VQE optimization loop (Python + Q#)
sweep_field.py           Python driver: sweeps the field to look for the phase transition
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
- **8 unit tests** (`@Test` + `Fact`): term-count correctness for both boundary conditions, the 2-qubit edge case, wrap-bond correctness, field-term correctness, and — most importantly — energy estimates matching the analytically-known exact value for both boundary conditions.

### `Ansatz.qs` — the parameterized trial-state circuit

This file implements the parameterized quantum circuit, or **ansatz**, used by the Variational Quantum Eigensolver. The ansatz prepares a trial quantum state whose parameters are adjusted by the classical optimizer until the measured energy is minimized. Because the quality of the ansatz directly affects the accuracy of VQE, this file is responsible for generating the quantum states that are evaluated throughout the optimization process.


- **`RealAmplitudesAnsatz(theta, qubits, reps)`** — a hardware-efficient ansatz: `reps` repeated layers, each layer applying one tunable `Ry` rotation per qubit, followed by a chain of `CNOT`s linking each qubit to the next. `theta` holds `reps * numQubits` angles total (one per qubit per layer). Marked `is Adj + Ctl` so it composes cleanly with other Q# operations that need adjoint/controlled variants.
- The name follows Qiskit's `RealAmplitudes`, but the two are not identical. Qiskit's version applies one more rotation layer after the final entangling block, so its parameter count is `(reps + 1) * numQubits` rather than `reps * numQubits`. The structure is otherwise the same, and we checked that the missing layer does not reduce accuracy at the system sizes used here.
- Using only `Ry` rotations is justified by the physics rather than chosen for convenience. The Transverse-Field Ising Hamiltonian is stoquastic, which means its ground state can always be written using real, non-negative amplitudes. A circuit built from `Ry` rotations produces exactly those states, so it is expressive enough for this problem and no parameters are spent on complex phases that the answer does not need.
- This is the "how deep is the circuit" knob for the project's comparison: varying `reps` directly varies circuit depth, gate count, and (expected) accuracy — the sweep the resource-estimation and results sections are built around. The depth required grows with system size. A 2-spin chain is exact at `reps = 1`, a 3-spin chain at `reps = 2`, and a 4-spin chain at `reps = 3`.


### `VQE.qs` — connecting ansatz to Hamiltonian

This file combines the ansatz and Hamiltonian into a single quantum routine that can be called by the classical optimizer. Rather than performing the optimization itself, it prepares a parameterized quantum state, measures its energy with respect to the Hamiltonian, and returns that energy as a single value. During VQE, this operation is called hundreds or even thousands of times as the optimizer searches for the lowest-energy state.


- **`VQEEnergy(theta, numQubits, J, h, reps, periodic, shotsPerTerm)`** — builds the Hamiltonian, builds the ansatz-prepared state from `theta`, and returns the estimated energy via `EstimateEnergy`. This is the single function the classical optimizer calls, over and over, with different `theta` each time.

### `exact.py` — classical ground-truth baseline

VQE is an approximate algorithm, so its results need to be compared against an exact solution whenever possible. This file provides that classical benchmark by constructing the same Hamiltonian as a NumPy matrix and solving it using exact diagonalization. The resulting ground-state energy is used to evaluate how accurately the quantum algorithm performs.


- **`to_matrix(H)`** — normalizes a Hamiltonian (however it's represented) to a plain NumPy array.
- **`exact_ground_energy(H)`** — lowest eigenvalue via `np.linalg.eigvalsh`.
- **`exact_ground_state(H)`** — lowest eigenvalue *and* its eigenvector.
- **`exact_spectrum(H)`** — every energy level, sorted from lowest to highest.
- **`tfim_open_chain_energy(n, J, h)`** — the ground energy of the open chain, computed without building the `2^n` matrix at all. The Transverse-Field Ising chain can be mapped onto a system of non-interacting fermions, and for an open chain the ground energy then works out to minus the sum of the singular values of a small `n` by `n` matrix. This agrees with full diagonalization to machine precision, and because the matrix stays small it remains fast for hundreds of spins. That makes it possible to show how the results behave at system sizes far beyond what the quantum simulation can reach.
- Includes a runnable self-check (`__main__` block) verifying two known cases by hand: `H = -X` (ground energy `-1`) and the 2-spin Ising Hamiltonian `H = -Z₀Z₁ - X₀ - X₁` (ground energy `-√5 ≈ -2.236`), which is the textbook sanity-check value for this model.

### `test_exact.py` — unit tests for the classical baseline

This file contains the Python unit tests for the classical solver. Since the exact solution is used as the benchmark throughout the project, it is important to verify that the classical implementation is producing the correct results before comparing it against the quantum algorithm.


6 pytest tests: single-spin `-X` and `Z` ground energies, the known `-√5` two-spin value, spectrum length/ordering (2 qubits → 4 levels, sorted ascending), the eigen-equation check (`H·ψ = E·ψ` for the returned ground state), and a trivial diagonal-matrix case as a basic correctness anchor.

### `run_vqe.py` — the full VQE optimization loop (Python host + Q#)

This file is the main driver for the project and implements the hybrid quantum-classical optimization loop. It connects the Python optimizer to the Q# quantum routines by repeatedly proposing new ansatz parameters, calling the quantum energy evaluation, and updating the parameters until the estimated energy converges. This script is responsible for executing the complete VQE algorithm from start to finish.


- **`load_qsharp()`** — compiles `Ansatz.qs`, `Hamiltonian.qs`, `VQE.qs` into the running Python session via `qdk.qsharp`.
- **`energy(theta, n, J, h, reps, periodic, shots)`** — thin wrapper that formats `theta` into a Q# call to `VQEEnergy` and returns the result.
- **`nft(objective, start, maxiter)`** — our implementation of the NFT optimizer, also known as Rotosolve. Because the ansatz contains only `Ry` gates, the energy is an exact sine wave as a function of any single angle when the others are held fixed. Three measurements are enough to determine that sine wave completely, so the best value of each angle can be solved for directly instead of being approached in small steps. This also means the optimizer has no learning rate or step size that needs tuning.
- **`one_run(...)`** — a single optimization attempt. Passing `method="NFT"` uses the routine above, and any other value is handed to `scipy.optimize.minimize`, with `COBYLA` used as the baseline for comparison. The best angles found are then re-measured using more shots. This is necessary because the lowest sampled value seen during optimization is statistically biased downward by noise, so reporting it directly would make the result look more accurate than it is.
- **`run_vqe(...)`** — wraps `one_run` with **multiple random restarts** (default 5), since either optimizer can settle in a local minimum depending on its starting angles; picks the best restart, then does one final high-shot remeasurement of the winning angles for the reported number.
- **`at(op, i, n)`** and **`tfim_matrix(n, J, h, periodic)`** — build the *same* Hamiltonian as `Hamiltonian.qs`, but as a plain NumPy matrix, so `exact_ground_energy` can be run on it for comparison. (Worth double-checking these two Hamiltonian constructions — the Q# one and this NumPy one — stay in sync if either changes.)
- **`__main__` block** — runs VQE for `n = 2, 3`, compares against the exact floor, and prints a table of VQE energy, exact energy, absolute error, and number of energy evaluations used. It also saves a convergence plot per system size.

### `estimate_resources.py` — resource estimation

This file analyzes the quantum resources required to run the VQE circuit using the Microsoft Quantum Development Kit Resource Estimator. Rather than evaluating the accuracy of the algorithm, it estimates how expensive the quantum circuit would be on a fault-tolerant quantum computer by reporting gate counts, qubit counts, runtime estimates, and other hardware metrics. These results are used to study how the cost of the algorithm scales as the problem size increases.


- **`load_qsharp(filename)`** — compiles a single `.qs` file so its operations are callable.
- **`estimate_resources(entry_expr, params=None)`** — calls Q#'s built-in Resource Estimator on a Q# call expression (e.g. `"VQEEnergy(...)"`) and pulls out: logical qubit count, rotation-gate count, T-gate count, measurement count, physical qubit count, runtime, code distance, and number of T-factories.
- **`print_estimate(...)`** — pretty-prints one estimate.
- **`__main__` block** — sweeps system size `n ∈ {2, 3, 4}` and ansatz depth `reps ∈ {1, 2, 3}`, running the resource estimator on `VQEEnergy` at each combination. This is exactly the sweep the project's "circuit depth vs. resource cost" comparison needs.

### `Observables.qs` — measuring the physical properties of the state

This file measures physical properties of the quantum state other than its energy. The energy tells us how well the optimization worked, but it does not describe what the spins are actually doing, so this file provides the quantities that do. It reuses the `EstimateEnergy` operation from `Hamiltonian.qs` by passing it a custom list of Pauli terms instead of the Hamiltonian itself, which means no new measurement machinery is needed.

- **`TransverseTerms(numQubits)`** — builds the term list for the average transverse magnetisation, one `X` term per spin, each weighted by `1/numQubits`.
- **`CouplingTerms(numQubits, periodic)`** — builds the term list for the average neighbour correlation, one `ZZ` term per bond, using the same boundary condition as the Hamiltonian.
- **`MeasureTransverse(...)`** and **`MeasureCoupling(...)`** — prepare the ansatz state from `theta` and return those two averages.

These two quantities are what identify the phase of the system. When the field is strong the transverse magnetisation is close to 1 and the neighbour correlation is close to 0, and when the field is weak the two values swap.

The longitudinal magnetisation is the order parameter normally used for the Ising model, but it is not measured here. On a finite chain its value is exactly zero at every field strength, because the ground state is an equal superposition of all spins up and all spins down and the two contributions cancel. Plotting it would produce a flat line at zero regardless of the physics.

### `collect_results.py` — running the sweep once and saving the numbers

This file runs the full sweep over system size and ansatz depth a single time and writes every result to `results.json`. Because each VQE run takes several minutes, re-running the algorithm every time a plot or a table needs adjusting is impractical. Saving the results once means the plots can be rebuilt from stored data in seconds, and it also means the numbers reported in the writeup stay fixed instead of changing slightly on every run.

For each configuration it stores the VQE energy, the exact energy, the error, the number of parameters, the full convergence history, the final angles, the elapsed time, and the counts from the resource estimator. The file is written after each configuration finishes, so a sweep that is interrupted still keeps everything completed up to that point.

### `sweep_field.py` — sweeping the field to find the phase transition

This file keeps the system size fixed and varies the field strength instead. The Ising model is known for having a phase transition as the field is changed. At each field strength the script optimizes the circuit and then measures the energy along with both observables.

Two choices in this script are worth explaining. The sweep runs downward from a strong field, and each solution is reused as the starting point for the next field value. This is done because the strong-field ground state is easy for the ansatz to represent, so starting there and moving gradually keeps the optimizer close to a good solution throughout. The sweep also uses the periodic ring rather than the open chain, because at four spins an open chain has a large finite-size effect that moves the transition signature far away from where it belongs.

### `requirements.txt`

This file lists all of the Python dependencies required to run the classical portion of the project. Installing these packages ensures that the optimization routines, testing framework, numerical computations, resource estimation scripts, and the Python–Q# interface all function correctly without additional manual setup.

Lists the Python dependencies needed to run the host-side scripts (the `qdk.qsharp` package, `numpy`, `scipy`, `pytest`, based on what `run_vqe.py`, `exact.py`, `estimate_resources.py`, and `test_exact.py` import). `qdk` is pinned because the resource estimator we call is on a deprecated code path that a future release may remove.

---

## Running it

```
py -3.13 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

```
.venv\Scripts\python.exe -m pytest -q            # unit tests for the classical baseline
.venv\Scripts\python.exe exact.py                # classical baseline self-check
.venv\Scripts\python.exe estimate_resources.py   # circuit cost across sizes and depths
.venv\Scripts\python.exe run_vqe.py              # full VQE run, takes a few minutes
.venv\Scripts\python.exe collect_results.py      # the whole sweep, saved to results.json
.venv\Scripts\python.exe sweep_field.py          # the field sweep, saved to sweep.json
```

---

## Results

### Accuracy against the classical answer

The VQE energy matches the exact ground-state energy at every system size tested. The total energy grows as spins are added, so the error per site is included as well, since that is the quantity that can be compared fairly across different chain lengths. The values below use NFT at two ansatz layers, and come from the saved sweep in `results.json`.

| Spins | VQE | Exact | Error | Error per site |
|-------|-----|-------|-------|----------------|
| 2 | -2.2349 | -2.2361 | 0.0012 | 0.0006 |
| 3 | -3.4960 | -3.4940 | 0.0021 | 0.0007 |
| 4 | -4.7513 | -4.7588 | 0.0075 | 0.0019 |

The 3-spin value sits slightly below the exact energy. The energy is estimated from a finite number of measurements, so the result is a random quantity centred on the true value, and a difference this small is well within the expected statistical spread. The limitations section below explains why this happens and why it does not contradict the variational principle.

### The classical optimizer mattered more than the circuit

Our first working version used COBYLA, and the error grew quickly as the chain got longer. The obvious explanation was that the ansatz was too shallow to represent the ground state at larger sizes, so we tested that directly by optimizing the same circuit with no sampling noise. The circuit reached the exact energy at all three sizes, which showed that the ansatz was not the limitation.

The actual cause was the optimizer. COBYLA works by building a simple model of the energy over a shrinking region, and once that region becomes small enough, the real differences in energy across it are smaller than the noise from finite sampling. At that point the model is fitting noise, and COBYLA treats this as convergence and stops early. This becomes worse as the number of parameters grows, which is why the problem appeared at larger system sizes.

Replacing COBYLA with NFT removed the problem. Since the ansatz uses only `Ry` gates, the energy is an exact sine wave in each angle, and NFT determines the best value of each angle in closed form from three measurements rather than searching for it. There is no small region for noise to interfere with.

We ran both optimizers at every combination of system size and ansatz depth, so the comparison below covers nine configurations. Times are for the whole run including restarts.

| Spins | Layers | NFT error | NFT time | COBYLA error | COBYLA time |
|-------|--------|-----------|----------|--------------|-------------|
| 2 | 1 | 0.0031 | 211 s | 0.0188 | 40 s |
| 2 | 2 | 0.0012 | 269 s | 0.0117 | 65 s |
| 2 | 3 | 0.0080 | 318 s | 0.0097 | 86 s |
| 3 | 1 | 0.0135 | 442 s | 0.0117 | 98 s |
| 3 | 2 | 0.0021 | 568 s | 0.0432 | 170 s |
| 3 | 3 | 0.0030 | 738 s | 0.0044 | 308 s |
| 4 | 1 | 0.0337 | 746 s | 0.9934 | 165 s |
| 4 | 2 | 0.0075 | 1099 s | 0.0765 | 475 s |
| 4 | 3 | 0.0092 | 1323 s | 0.8888 | 734 s |

NFT is more accurate in eight of the nine configurations. The one exception is the shallowest 3-spin circuit, where the two are effectively tied.

In two of the nine runs, both at four spins, COBYLA finished with an error close to 1.0, meaning it failed to find the ground state. NFT never failed this way, and its worst result across all nine runs was an error of 0.034. An optimizer that occasionally returns a completely wrong answer is harder to work with than one that is consistently a little imprecise, because there is no way to tell from the output alone which kind of run you got.

This accuracy comes at a cost in time. NFT took between two and five times longer than COBYLA, with a median of 568 seconds against 165 seconds. This is expected, since NFT spends three measurements on every parameter in every sweep, while COBYLA stops as soon as its own convergence test is satisfied. Part of the reason COBYLA is faster is that it stops early, which is also the reason it sometimes stops at the wrong answer.

### Finding the phase transition

Sweeping the field from 2.0 down to 0 on a 4-spin ring produces the change of phase directly. The transverse magnetisation falls from 0.919 to 0.000 as the field is removed, and the neighbour correlation rises from 0.294 to 1.000 over the same range. At strong field each spin follows the field independently and neighbouring spins are almost uncorrelated. At weak field the coupling wins and every neighbouring pair is locked together.

The point where the two curves cross is what identifies the transition, and this is where the choice of a ring rather than an open chain matters. The measured crossing is at `h/J = 1.004`, against a true critical point of 1.000.

At `h = J` the periodic Ising chain is self-dual, meaning the model maps onto itself with the coupling and the field exchanged. Because that mapping also exchanges the two quantities we are measuring, they are forced to be equal at exactly `h = J`, for any number of spins. Exact diagonalization confirms this: the crossing sits at 1.0000 for rings of 3, 4, 5, and 6 spins with no drift at all. On an open chain the same crossing is at 0.649, 0.729, 0.780, and 0.815 for those sizes, approaching 1 only slowly, so a 4-spin open chain would have placed the transition around 0.73 and required an explanation for the discrepancy.

The energy error across the sweep averages 0.0115 and is largest near the transition, which is the region where the optimization is hardest.

### Resource cost on real hardware

The resource estimator reports what would be required to run the circuit on a fault-tolerant quantum computer. The difference between the logical and physical qubit counts is the overhead of quantum error correction.

| Spins | Logical qubits | Rotation gates | Physical qubits |
|-------|----------------|----------------|-----------------|
| 2 | 3 | 12 | 225,544 |
| 3 | 4 | 30 | 323,070 |
| 4 | 5 | 56 | 391,556 |

These numbers illustrate the scale of that overhead. The estimator assumes a fault-tolerant machine, whereas VQE is an algorithm designed for near-term hardware that does not have error correction. The T-gate count also depends on a rotation-synthesis precision that we selected, so it is not a fixed property of the circuit.

---

## Limitations

- Energies are estimated from a finite number of measurements, so every evaluation carries some statistical noise.
- The variational principle guarantees that the true expectation value of the energy cannot be lower than the ground-state energy. A sampled estimate, however, is a random quantity and can fall below it by chance. For this reason we re-measure the final angles using a large number of shots rather than reporting the lowest value observed during optimization.
- Optimization is most difficult near the phase transition. The error there is a few times larger than at either end of the field sweep.
- The reference paper uses a periodic ring of twelve spins at a different field strength, so our energies cannot be compared directly against the values reported there.
- Exact diagonalization is used only as a benchmark and is not part of the algorithm. It is limited to small systems, since the matrix it works with doubles in size with every spin added.


