from pathlib import Path

import qdk.qsharp as qs


def load_qsharp(filename):
    # compile a .qs file from this folder so we can estimate its operations
    path = Path(__file__).parent / filename
    qs.eval(path.read_text(encoding="utf-8"))


def estimate_resources(entry_expr, params=None):
    # ask the estimator about a Q# call like "MiniAnsatz()"
    result = qs.estimate(entry_expr, params) if params else qs.estimate(entry_expr)
    logical = result["logicalCounts"]
    physical = result["physicalCounts"]
    pretty = result["physicalCountsFormatted"]
    return {
        "logical_qubits": logical["numQubits"],
        "rotation_count": logical["rotationCount"],
        "t_count": logical["tCount"],
        "measurement_count": logical["measurementCount"],
        "physical_qubits": physical["physicalQubits"],
        "physical_qubits_pretty": pretty["physicalQubits"],
        "runtime_pretty": pretty["runtime"],
        "code_distance": result["logicalQubit"]["codeDistance"],
        "num_tfactories": physical["breakdown"]["numTfactories"],
    }


def print_estimate(entry_expr, params=None):
    r = estimate_resources(entry_expr, params)
    print(f"Estimate for {entry_expr}")
    print(f"  logical qubits   {r['logical_qubits']}")
    print(f"  rotation gates   {r['rotation_count']}")
    print(f"  T gates          {r['t_count']}")
    print(f"  measurements     {r['measurement_count']}")
    print(f"  physical qubits  {r['physical_qubits']:,} ({r['physical_qubits_pretty']})")
    print(f"  runtime          {r['runtime_pretty']}")
    print(f"  code distance    {r['code_distance']}")
    print(f"  T-factories      {r['num_tfactories']}")
    return r


if __name__ == "__main__":
    for name in ("Ansatz.qs", "Hamiltonian.qs", "VQE.qs"):
        load_qsharp(name)

    # measure the real VQE circuit at a few system sizes
    for n in (2, 3, 4):
        theta = ", ".join(["0.5"] * (2 * n))
        print_estimate(f"VQEEnergy([{theta}], {n}, 1.0, 1.0, 1)")
        print()
