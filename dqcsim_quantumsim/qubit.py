from dqcsim.plugin import *
from dqcsim.common import *

class Qubit:
    def __init__(self, qsi, qubit_ref):
        super().__init__()
        self.qsi = qsi

        # Upstream/DQCsim qubit reference.
        self.qubit_ref = qubit_ref

        # QuantumSim SparseDM qubit index associated with this qubit. This
        # is None between measurements/allocs and the first gate.
        self._qs_ref = None

        # If this qubit is not in the SDM, remember the most recent
        # measurement for when we have to add it again. None is used for
        # undefined.
        self.classical = 0

    @property
    def qs_ref(self):
        return self._qs_ref

    @qs_ref.setter
    def qs_ref(self, new_qs_ref):
        if self._qs_ref == new_qs_ref:
            return
        if self._qs_ref is not None:
            self.qsi.free_qs_qubits.add(self._qs_ref)
            self.qsi.live_qs_qubits.remove(self._qs_ref)
            self._qs_ref = None
        if new_qs_ref is not None:
            self.qsi.free_qs_qubits.remove(new_qs_ref)
            self.qsi.live_qs_qubits.add(new_qs_ref)
            self._qs_ref = new_qs_ref

    def measure(self):
        """Measure this qubit in the Z basis."""

        # If this qubit is live within the SDM, observe it.
        if self.qs_ref is not None:

            # Get the measurement probabilities for this qubit.
            p0, p1 = self.qsi.sdm.peak_measurement(self.qs_ref)
            trace = p0 + p1

            # This is the total probability for the event up to this point,
            # including all past measurements, so p0 and p1 might add up
            # to less than one.
            p1 /= p0 + p1
            self.classical = int(self.qsi.random_float() < p1)

            # Project the measurement.
            self.qsi.sdm.project_measurement(self.qs_ref, int(bool(self.classical)))

            # Renormalize when the trace becomes too low to prevent numerical
            # problems (we don't use the trace for anything in this plugin).
            if trace < 1e-10:
                self.qsi.debug('renormalizing state density matrix, trace was {}...', trace)
                self.qsi.sdm.renormalize()

            # The qubit is now no longer relevant in the SDM, at least until
            # the next gate is applied. So we can take it out.
            self.qs_ref = None

        return Measurement(self.qubit_ref, self.classical)

    def prep(self):
        """Put this qubit in the |0> state."""

        # Measure ourself to make ourselves classical.
        self.measure()

        # Set our state to 0.
        self.classical = 0

    def ensure_in_sdm(self):
        """Make sure this qubit is represented in the SDM. Opposite of
        measure(), in a way. This must be called before a gate can be
        applied to the qubit."""
        if self.qs_ref is None:

            # Find a free SDM index.
            try:
                qs_ref = next(iter(self.qsi.free_qs_qubits))
            except StopIteration:
                raise RuntimeError(
                    'Too many qubits in use! Max is currently fixed to {}'
                    .format(self.qsi.MAX_QUBITS))

            # Claim the index.
            self.qs_ref = qs_ref

            # Make sure the SDM has the right bit value set.
            assert self.qs_ref in self.qsi.sdm.classical
            self.qsi.sdm.classical[self.qs_ref] = int(bool(self.classical))
