namespace VQE.Ansatz {

    open Microsoft.Quantum.Intrinsic;

    operation RealAmplitudesAnsatz(
        theta : Double[],
        qubits : Qubit[],
        reps : Int
    ) : Unit is Adj + Ctl {

        let n = Length(qubits);

        for rep in 0 .. reps - 1 {

            // rotation layer, one tunable angle per qubit
            for i in 0 .. n - 1 {
                Ry(theta[rep * n + i], qubits[i]);
            }

            // entanglement layer, links each qubit to the next
            for i in 0 .. n - 2 {
                CNOT(qubits[i], qubits[i + 1]);
            }
        }
    }
}
