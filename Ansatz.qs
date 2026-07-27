namespace VQE.Ansatz {

    open Microsoft.Quantum.Intrinsic;

    operation RealAmplitudesAnsatz(
        theta : Double[],
        qubits : Qubit[]
    ) : Unit is Adj + Ctl {

        let n = Length(qubits);
       
        // First Rotation Layer
        for i in 0 .. n - 1 {
            Ry(theta[i], qubits[i]);
        }

        // Linear Entanglement Layer       
        for i in 0 .. n - 2 {
            CNOT(qubits[i], qubits[i + 1]);
        }
        
        // Second Rotation Layer
        let offset = n;

        for i in 0 .. n - 1 {
            Ry(theta[offset + i], qubits[i]);
        }
        
        // Second Entanglement Layer
        for i in 0 .. n - 2 {
            CNOT(qubits[i], qubits[i + 1]);
        }
    }
}