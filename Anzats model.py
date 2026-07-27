from qiskit import QuantumCircuit


def build_real_amplitudes_ansatz(num_qubits, theta):
    """
   this is kinda 
    """

    qc = QuantumCircuit(num_qubits)

    ###################################################
    # First Rotation Layer
    ###################################################

    for i in range(num_qubits):
        qc.ry(theta[i], i)

    ###################################################
    # First Entanglement Layer
    ###################################################

    for i in range(num_qubits - 1):
        qc.cx(i, i + 1)

    ###################################################
    # Second Rotation Layer
    ###################################################

    offset = num_qubits

    for i in range(num_qubits):
        qc.ry(theta[offset + i], i)

    ###################################################
    # Second Entanglement Layer
    ###################################################

    for i in range(num_qubits - 1):
        qc.cx(i, i + 1)

    return qc
