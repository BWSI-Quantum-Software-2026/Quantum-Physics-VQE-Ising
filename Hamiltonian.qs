// Hamiltonian.qs

import Std.Diagnostics.*;
import Std.Convert.*;
import Std.Math.*;

struct HamiltonianTerm {
    Coefficient : Double,
    Paulis : Pauli[],
    Targets : Int[]
}



function GenerateTFIMHamiltonian(numQubits : Int, J : Double, h : Double) : HamiltonianTerm[] {
    Fact(numQubits >= 1, "numQubits must be at least 1.");

    mutable terms : HamiltonianTerm[] = [];

    for i in 0..numQubits - 2 {
        set terms += [HamiltonianTerm(-J, [PauliZ, PauliZ], [i, i + 1])];
    }


    for i in 0..numQubits - 1 {
        set terms += [HamiltonianTerm(-h, [PauliX], [i])];
    }

    return terms;
}


function SelectTargets(qubits : Qubit[], targets : Int[]) : Qubit[] {
    mutable selected : Qubit[] = [];
    for index in targets {
        set selected += [qubits[index]];
    }
    return selected;
}


operation MeasureTermEigenvalue(qubits : Qubit[], term : HamiltonianTerm) : Double {
    let targetQubits = SelectTargets(qubits, term.Targets);
    let result = Measure(term.Paulis, targetQubits);
    return result == Zero ? 1.0 | -1.0;
}


operation EstimateEnergy(
    prepareState : (Qubit[] => Unit),
    numQubits : Int,
    terms : HamiltonianTerm[],
    shotsPerTerm : Int
) : Double {
    mutable energy = 0.0;

    for term in terms {
        mutable eigenvalueSum = 0.0;

        for _ in 1..shotsPerTerm {
            use qubits = Qubit[numQubits];
            prepareState(qubits);
            set eigenvalueSum += MeasureTermEigenvalue(qubits, term);
            ResetAll(qubits);
        }

        let averageEigenvalue = eigenvalueSum / IntAsDouble(shotsPerTerm);
        set energy += term.Coefficient * averageEigenvalue;
    }

    return energy;
}


operation PrepareAllZero(qubits : Qubit[]) : Unit {
    // No-op: qubits are already |0...0> after allocation.
}


operation VerifyHamiltonian() : Unit {
    let numQubits = 4;
    let J = 1.0;
    let h = 0.5;
    let terms = GenerateTFIMHamiltonian(numQubits, J, h);

    Fact(Length(terms) == (numQubits - 1) + numQubits,
        $"Expected {(numQubits - 1) + numQubits} terms, got {Length(terms)}.");


    let expectedEnergy = -J * IntAsDouble(numQubits - 1);


    let estimatedEnergy = EstimateEnergy(PrepareAllZero, numQubits, terms, 2000);

    Message($"Expected energy (analytic): {expectedEnergy}");
    Message($"Estimated energy (sampled): {estimatedEnergy}");

    Fact(AbsD(estimatedEnergy - expectedEnergy) < 0.2,
        "Estimated energy deviates from the analytic value by more than the expected sampling noise.");

    Message("VerifyHamiltonian passed.");
}

@EntryPoint()
operation Main() : Unit {
    VerifyHamiltonian();
}
