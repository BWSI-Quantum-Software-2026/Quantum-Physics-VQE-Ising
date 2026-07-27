namespace VQE.Ansatz {

    open Microsoft.Quantum.Intrinsic;

    operation RealAmplitudesAnsatz(
        theta : Double[],
        qubits : Qubit[],
        reps : Int
    ) : Unit is Adj + Ctl {

        let n = Length(qubits);
        mutable thetaIndex = 0;
        
        for rep in 1 .. reps{ 
             // Rotation layer
            for i in 0 .. n - 1 {

                Ry(theta[thetaIndex], qubits[i]);

                set thetaIndex += 1;
            }

            // Linear entanglement
            for i in 0 .. n - 2 {

                CNOT(qubits[i], qubits[i + 1]);

            }
        
        }
    }
}
