from dqcsim_quantumsim.qubit import Qubit
from dqcsim.plugin import *
from dqcsim.common import *

@plugin("QuantumSim interface", "Jeroen van Straten", "0.0.1")
class QuantumSimInterface(Backend):

    # QuantumSim's SparseDM object doesn't support adding or removing qubits.
    # However, any qubits that haven't been entangled yet are purely classical.
    # Therefore, we can just allocate a large number of qubits at startup and
    # use those when we need them. This is that large number.
    MAX_QUBITS = 1000

    def __init__(self):
        super().__init__()

        # quantumsim.ptm module reference.
        self.ptm = None

        # quantumsim.sparsedm.SparseDM object representing the system.
        self.sdm = None

        # numpy module reference.
        self.np = None

        # Indices of qubits that are free/live within self.sdm.
        self.free_qs_qubits = set(range(self.MAX_QUBITS))
        self.live_qs_qubits = set()

        # Qubit data for each upstream qubit.
        self.qubits = {}

    def handle_init(self, *_a, **_k):
        # Loading QuantumSim can take some time, so defer to initialize
        # callback. We also have logging at that point in time, so it should
        # provide a nice UX.
        self.debug('Trying to load QuantumSim...')
        import quantumsim.ptm as ptm
        self.ptm = ptm
        import quantumsim.sparsedm as sdm
        self.sdm = sdm.SparseDM(self.MAX_QUBITS)
        import numpy as np
        self.np = np
        self.info('QuantumSim loaded {}using CUDA acceleration', '' if sdm.using_gpu else '*without* ')

    def handle_allocate(self, qubit_refs, *_a, **_k):
        for qubit_ref in qubit_refs:
            self.qubits[qubit_ref] = Qubit(self, qubit_ref)

    def handle_free(self, qubit_refs, *_a, **_k):
        for qubit_ref in qubit_refs:
            qubit = self.qubits.pop(qubit_ref)

            # Measure the qubit to make sure it's freed in the SDM.
            qubit.measure()

    def handle_measurement_gate(self, qubit_refs, *_a, **_k):
        measurements = []
        for qubit_ref in qubit_refs:
            qubit = self.qubits[qubit_ref]
            measurements.append(qubit.measure())
        return measurements

    def handle_unitary_gate(self, qubit_refs, pauli_matrix, *_a, **_k):
        if len(qubit_refs) == 1:

            # Single-qubit gate. Unpack the qubit from the list.
            qubit_ref, = qubit_refs

            # Make sure the qubit is present in the SDM.
            qubit = self.qubits[qubit_ref]
            qubit.ensure_in_sdm()

            # Convert the incoming matrix to a numpy array.
            pauli_matrix = self.np.array([
                pauli_matrix[0:2],
                pauli_matrix[2:4]])

            # Convert the Pauli matrix to the corresponding ptm.
            ptm = self.ptm.single_kraus_to_ptm(pauli_matrix)

            # Apply the ptm.
            self.sdm.apply_ptm(qubit.qs_ref, ptm)

        elif len(qubit_refs) == 2:

            # Two-qubit gate. Unpack the qubits from the list.
            qubit_ref_a, qubit_ref_b = qubit_refs

            # Make sure the qubits are present in the SDM.
            qubit_a = self.qubits[qubit_ref_a]
            qubit_a.ensure_in_sdm()

            qubit_b = self.qubits[qubit_ref_b]
            qubit_b.ensure_in_sdm()

            # Convert the incoming matrix to a numpy array.
            pauli_matrix = self.np.array([
                pauli_matrix[0:4],
                pauli_matrix[4:8],
                pauli_matrix[8:12],
                pauli_matrix[12:16]])

            # Convert the Pauli matrix to the corresponding ptm.
            two_ptm = self.ptm.double_kraus_to_ptm(pauli_matrix)

            # Apply the ptm.
            self.sdm.apply_two_ptm(qubit_a.qs_ref, qubit_b.qs_ref, two_ptm)

        else:
            raise RuntimeError(
                'QuantumSim can only handle one- and two-qubit gates. ' +
                '{} is too many.'.format(len(qubit_refs)))

