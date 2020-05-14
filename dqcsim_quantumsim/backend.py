from dqcsim_quantumsim.qubit import Qubit
from dqcsim.plugin import *
from dqcsim.common import *

@plugin("QuantumSim interface", "Jeroen van Straten", "0.0.5")
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

    def handle_init(self, cmds):

        # Interpret commands.
        self.t1 = None
        self.t2 = None
        for cmd in cmds:
            if cmd.iface == 'quantumsim':
                if cmd.oper == 'error':
                    if 't1' in cmd:
                        self.t1 = cmd['t1']
                    if 't2' in cmd:
                        self.t2 = cmd['t2']
                else:
                    raise ValueError('Unknown ArbCmd %s.%s' % (cmd.iface, cmd.oper))

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

    def handle_allocate(self, qubit_refs, cmds):

        # Allow t1 and t2 to be overwritten on a per-qubit basis.
        t1 = self.t1
        t2 = self.t2
        for cmd in cmds:
            if cmd.iface == 'quantumsim':
                if cmd.oper == 'error':
                    if 't1' in cmd:
                        t1 = cmd['t1']
                    if 't2' in cmd:
                        t2 = cmd['t2']
                else:
                    raise ValueError('Unknown ArbCmd %s.%s' % (cmd.iface, cmd.oper))

        for qubit_ref in qubit_refs:
            self.qubits[qubit_ref] = Qubit(self, qubit_ref, t1=t1, t2=t2)

    def handle_free(self, qubit_refs):
        for qubit_ref in qubit_refs:
            qubit = self.qubits.pop(qubit_ref)

            # Measure the qubit to make sure it's freed in the SDM.
            qubit.measure()

    def handle_measurement_gate(self, qubit_refs, basis, arb):

        # Determine the projection method.
        if 'method' in arb:
            methods = arb['method']
        else:
            methods = 'random'
        if isinstance(methods, list):
            if len(methods) != len(qubit_refs):
                raise ValueError('method key does not have the right list size')
        elif isinstance(methods, str):
            methods = [methods] * len(qubit_refs)
        elif isinstance(methods, int):
            methods = [(methods >> i) & 1 for i in reversed(range(len(qubit_refs)))]
        else:
            raise ValueError('failed to parse method key')

        # Determine the hermitian of the basis for rotating back after
        # measuring.
        basis_hermetian = [
            basis[0].real - basis[0].imag * 1j,
            basis[2].real - basis[2].imag * 1j,
            basis[1].real - basis[1].imag * 1j,
            basis[3].real - basis[3].imag * 1j]

        # Perform the measurements.
        measurements = []
        for qubit_ref, method in zip(qubit_refs, methods):
            self.handle_unitary_gate([qubit_ref], basis_hermetian, None)
            qubit = self.qubits[qubit_ref]
            measurements.append(qubit.measure(method))
            self.handle_unitary_gate([qubit_ref], basis, None)
        return measurements

    def handle_prepare_gate(self, qubit_refs, basis, _arb):
        measurements = []
        for qubit_ref in qubit_refs:
            self.qubits[qubit_ref].prep()
            self.handle_unitary_gate([qubit_ref], basis, None)
        return measurements

    def handle_unitary_gate(self, qubit_refs, unitary_matrix, _arb):

        # Convert the incoming matrix to a numpy array.
        unitary_matrix = self.np.reshape(self.np.array(unitary_matrix), (2**len(qubit_refs),) * 2)

        # Print what we're doing.
        self.debug('gate on %s:\n%s' % (', '.join(map('q{}'.format, qubit_refs)), unitary_matrix))

        # Apply pending idling gates and get QuantumSim references.
        qs_refs = []
        for qubit_ref in qubit_refs:
            qubit = self.qubits[qubit_ref]
            qubit.ensure_in_sdm()
            qubit.apply_pending_error()
            qs_refs.append(qubit.qs_ref)

        # Convert the Z-basis unitary matrix to the corresponding Pauli transfer
        # matrix in QuantumSim's 0xy1 basis.
        tensor = self.ptm.single_tensor
        for _ in range(len(qubit_refs) - 1):
            tensor = self.np.kron(tensor, self.ptm.single_tensor)
        ptm = self.np.einsum("xab, bc, ycd, ad -> xy", tensor, unitary_matrix, tensor, unitary_matrix.conj()).real

        # Apply the ptm.
        if len(qubit_refs) == 1:
            self.sdm.apply_ptm(qs_refs[0], ptm)
        elif len(qubit_refs) == 2:
            self.sdm.apply_two_ptm(qs_refs[1], qs_refs[0], ptm)
        else:
            raise RuntimeError(
                'QuantumSim can only handle one- and two-qubit gates. ' +
                '{} is too many.'.format(len(qubit_refs)))

    def handle_advance(self, cycles):
        for qubit in self.qubits.values():
            qubit.idle(cycles)
