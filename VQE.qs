// # Sample
// Getting started
//
// # Description
// This is a minimal Q# program that can be used to start writing Q# code.

operation VQE() : Unit {

    mutable theta = [0.5, 0.7, 1.2, 0.3];

    let stepSize = 0.1;
    let maxIterations = 100;

    mutable bestEnergy = 10000000000.0;
    mutable bestTheta = theta;

    use qubits = Qubit[2];

    for i in 0 .. maxIterations - 1 {

        Ansatz(theta, qubits);
        let energy = MeasureEnergy(qubits);

        if energy < bestEnergy {
            set bestEnergy = energy;
            set bestTheta = theta;
        }

        let gradient = ComputeGradient(theta);

        set theta =
        UpdateParameters(
            theta,
            gradient,
            learningRate
        );

        ResetAll(qubits);

    }

    Message($"Best Energy = {bestEnergy}");


}
