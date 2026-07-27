// builds the ansatz state and returns its energy.
// this file only measures energy.

operation VQEEnergy(
    theta : Double[],
    numQubits : Int,
    J : Double,
    h : Double,
    reps : Int,
    periodic : Bool,
    shotsPerTerm : Int
) : Double {
    let terms = GenerateTFIMHamiltonian(numQubits, J, h, periodic);
    return EstimateEnergy(
        VQE.Ansatz.RealAmplitudesAnsatz(theta, _, reps), numQubits, terms, shotsPerTerm);
}
