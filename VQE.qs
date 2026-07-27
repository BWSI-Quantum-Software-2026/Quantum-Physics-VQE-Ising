// builds the ansatz state and returns its energy.
// this file only measures energy.

operation VQEEnergy(
    theta : Double[],
    numQubits : Int,
    J : Double,
    h : Double,
    shotsPerTerm : Int
) : Double {
    let terms = GenerateTFIMHamiltonian(numQubits, J, h);
    return EstimateEnergy(
        VQE.Ansatz.RealAmplitudesAnsatz(theta, _), numQubits, terms, shotsPerTerm);
}
