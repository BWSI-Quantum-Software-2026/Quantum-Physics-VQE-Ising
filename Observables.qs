// measures the two physics observables we plot across the field sweep

import Std.Convert.*;

function TransverseTerms(numQubits : Int) : HamiltonianTerm[] {
    // average transverse magnetisation, one X per spin
    mutable terms : HamiltonianTerm[] = [];
    let weight = 1.0 / IntAsDouble(numQubits);
    for i in 0..numQubits - 1 {
        set terms += [HamiltonianTerm(weight, [PauliX], [i])];
    }
    return terms;
}


function CouplingTerms(numQubits : Int, periodic : Bool) : HamiltonianTerm[] {
    // average neighbour correlation, one ZZ per bond
    mutable pairs : (Int, Int)[] = [];
    for i in 0..numQubits - 2 {
        set pairs += [(i, i + 1)];
    }
    if periodic and numQubits > 2 {
        set pairs += [(numQubits - 1, 0)];
    }

    let weight = 1.0 / IntAsDouble(Length(pairs));
    mutable terms : HamiltonianTerm[] = [];
    for (a, b) in pairs {
        set terms += [HamiltonianTerm(weight, [PauliZ, PauliZ], [a, b])];
    }
    return terms;
}


operation MeasureTransverse(
    theta : Double[],
    numQubits : Int,
    reps : Int,
    shotsPerTerm : Int
) : Double {
    return EstimateEnergy(
        VQE.Ansatz.RealAmplitudesAnsatz(theta, _, reps),
        numQubits, TransverseTerms(numQubits), shotsPerTerm);
}


operation MeasureCoupling(
    theta : Double[],
    numQubits : Int,
    reps : Int,
    periodic : Bool,
    shotsPerTerm : Int
) : Double {
    return EstimateEnergy(
        VQE.Ansatz.RealAmplitudesAnsatz(theta, _, reps),
        numQubits, CouplingTerms(numQubits, periodic), shotsPerTerm);
}
