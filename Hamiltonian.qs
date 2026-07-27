// Hamiltonian.qs


import Std.Diagnostics.*;
import Std.Convert.*;
import Std.Math.*;


struct HamiltonianTerm {
    Coefficient : Double,
    Paulis : Pauli[],
    Targets : Int[]
}


function GenerateTFIMHamiltonian(numQubits : Int, J : Double, h : Double, periodic : Bool) : HamiltonianTerm[] {
    Fact(numQubits >= 1, "numQubits must be at least 1.");

    mutable terms : HamiltonianTerm[] = [];


    for i in 0..numQubits - 2 {
        set terms += [HamiltonianTerm(-J, [PauliZ, PauliZ], [i, i + 1])];
    }


    if periodic and numQubits > 2 {
        set terms += [HamiltonianTerm(-J, [PauliZ, PauliZ], [numQubits - 1, 0])];
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


@EntryPoint()
operation Main() : Unit {
    let numQubits = 4;
    let J = 1.0;
    let h = 0.5;
    let shots = 2000;

    let openTerms = GenerateTFIMHamiltonian(numQubits, J, h, false);
    let ringTerms = GenerateTFIMHamiltonian(numQubits, J, h, true);

    let openEnergy = EstimateEnergy(PrepareAllZero, numQubits, openTerms, shots);
    let ringEnergy = EstimateEnergy(PrepareAllZero, numQubits, ringTerms, shots);

    Message($"Open chain   ({Length(openTerms)} terms): estimated <H> = {openEnergy}, analytic = {-J * IntAsDouble(numQubits - 1)}");
    Message($"Periodic ring ({Length(ringTerms)} terms): estimated <H> = {ringEnergy}, analytic = {-J * IntAsDouble(numQubits)}");
    Message($"Difference from the extra wrap-around bond: analytic gap = {-J}, observed gap = {ringEnergy - openEnergy}");
}


@Test()
function OpenChainTermCountIsCorrect() : Unit {
    let numQubits = 4;
    let terms = GenerateTFIMHamiltonian(numQubits, 1.0, 0.5, false);
    let expectedCount = (numQubits - 1) + numQubits;
    Fact(Length(terms) == expectedCount,
        $"Open chain: expected {expectedCount} terms, got {Length(terms)}.");
}


@Test()
function PeriodicRingTermCountIsCorrect() : Unit {
    let numQubits = 4;
    let terms = GenerateTFIMHamiltonian(numQubits, 1.0, 0.5, true);
    let expectedCount = numQubits + numQubits;
    Fact(Length(terms) == expectedCount,
        $"Periodic ring: expected {expectedCount} terms, got {Length(terms)}.");
}


@Test()
function TwoQubitPeriodicDoesNotDuplicateBond() : Unit {
    let openTerms = GenerateTFIMHamiltonian(2, 1.0, 0.5, false);
    let periodicTerms = GenerateTFIMHamiltonian(2, 1.0, 0.5, true);
    Fact(Length(openTerms) == Length(periodicTerms),
        $"Expected periodic=true to leave the term count unchanged for a 2-qubit chain (both {Length(openTerms)}), got {Length(periodicTerms)} for periodic.");
}


@Test()
function WrapBondConnectsLastQubitToFirst() : Unit {
    let numQubits = 5;
    let terms = GenerateTFIMHamiltonian(numQubits, 1.0, 0.5, true);

    
    let wrapTerm = terms[numQubits - 1];

    Fact(wrapTerm.Targets[0] == numQubits - 1 and wrapTerm.Targets[1] == 0,
        $"Expected wrap bond to connect qubit {numQubits - 1} back to qubit 0, got {wrapTerm.Targets[0]} and {wrapTerm.Targets[1]}.");
    Fact(wrapTerm.Paulis[0] == PauliZ and wrapTerm.Paulis[1] == PauliZ,
        "Wrap bond should use PauliZ on both qubits, like every other ZZ bond.");
}


@Test()
function FieldTermUsesCorrectCoefficientAndBasis() : Unit {
    let J = 1.0;
    let h = 0.75;
    let numQubits = 3;
    let terms = GenerateTFIMHamiltonian(numQubits, J, h, false);

    
    let firstFieldTerm = terms[numQubits - 1];

    Fact(AbsD(firstFieldTerm.Coefficient - (-h)) < 1e-9,
        $"Expected field coefficient {-h}, got {firstFieldTerm.Coefficient}.");
    Fact(Length(firstFieldTerm.Paulis) == 1 and firstFieldTerm.Paulis[0] == PauliX,
        "Field term should act with a single PauliX.");
}


@Test()
operation OpenChainEnergyMatchesAnalyticValueForAllZeroState() : Unit {
    let numQubits = 4;
    let J = 1.0;
    let h = 0.5;
    let terms = GenerateTFIMHamiltonian(numQubits, J, h, false);
    let expectedEnergy = -J * IntAsDouble(numQubits - 1);
    let estimatedEnergy = EstimateEnergy(PrepareAllZero, numQubits, terms, 2000);

    Fact(AbsD(estimatedEnergy - expectedEnergy) < 0.2,
        $"Open chain |0...0> energy estimate {estimatedEnergy} deviates from analytic value {expectedEnergy} by more than expected sampling noise.");
}


@Test()
operation PeriodicRingEnergyMatchesAnalyticValueForAllZeroState() : Unit {
    let numQubits = 4;
    let J = 1.0;
    let h = 0.5;
    let terms = GenerateTFIMHamiltonian(numQubits, J, h, true);
    let expectedEnergy = -J * IntAsDouble(numQubits);
    let estimatedEnergy = EstimateEnergy(PrepareAllZero, numQubits, terms, 2000);

    Fact(AbsD(estimatedEnergy - expectedEnergy) < 0.2,
        $"Periodic ring |0...0> energy estimate {estimatedEnergy} deviates from analytic value {expectedEnergy} by more than expected sampling noise.");
}



@Test()
operation PeriodicHamiltonianHasCorrectNumberOfTerms() : Unit {
    let terms = GenerateTFIMHamiltonian(4, 1.0, 0.5, true);

    Message($"Number of terms = {Length(terms)}");

    Fact(
        Length(terms) == 8,
        "Expected 8 Hamiltonian terms."
    );
}
