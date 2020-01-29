from dqcsim.plugin import *
from dqcsim.common import *
import struct

class Qubit:
    def __init__(self, qsi, qubit_ref, t1=None, t2=None):
        """Construct a qubit wrapper for the given `QuantumSimInterface`,
        upstream qubit reference, and t1/t2 times.

        t1 is the amplitude damping time as measured in free decay experiments,
        t2 is the phase damping time as measured in ramsey or echo experiments.
        If unspecified or `None`, they are made infinite. Their unit is DQCsim
        cycles."""
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

        # Calculate error parameters.
        if t1 is None:
            t1 = qsi.np.inf
        if t2 is None:
            t2 = qsi.np.inf
        if qsi.np.allclose(t2, 2 * t1):
            t_phi = qsi.np.inf
        else:
            t_phi = 1 / (1 / t2 - 1 / (2 * t1)) / 2
        self.t1 = t1
        self.t2 = t2
        self.t_phi = t_phi

        # Number of cycles this qubit has been idle.
        self._idle_time = 0

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

    def measure(self, method='random'):
        """Measure this qubit in the Z basis.

        The following methods are supported:

         - 'random': calculate the probabilities for measuring a 0 vs a 1 and
           use DQCsim's PRNG to select one.
         - 'probable': calculate the probabilities for measuring a 0 vs a 1 and
           select the most probable.
         - 1: project the measurement to one regardless of probability.
         - 0: project the measurement to zero regardless of probability.

        The latter two return an error if the probability for projecting to the
        requested state is zero."""

        # Apply any pending error before measuring.
        self.apply_pending_error()

        # If this qubit is live within the SDM, observe it.
        if self.qs_ref is not None:

            # Get the measurement probabilities for this qubit.
            p0, p1 = self.qsi.sdm.peak_measurement(self.qs_ref)
            trace = p0 + p1
            p0 /= trace
            p1 /= trace

            # Select the value we're going to project to.
            if method == 'random':
                self.classical = int(self.qsi.random_float() < p1)
            elif method == 'probable':
                self.classical = int(p1 > p0)
            elif method == 0:
                if p0 < 1.e-20:
                    raise ValueError('cannot project to 0, probability too low at %f' % p0)
                self.classical = 0
            elif method == 1:
                if p1 < 1.e-20:
                    raise ValueError('cannot project to 1, probability too low at %f' % p1)
                self.classical = 1
            else:
                raise ValueError('unknown measurement method %r' % method)

            # Project the measurement.
            self.qsi.sdm.project_measurement(self.qs_ref, int(bool(self.classical)))

            # Determine and process the probability of this measurement.
            if self.classical:
                p = p1
            else:
                p = p0
            trace *= p

            # Renormalize when the trace becomes too low to prevent numerical
            # problems (we don't use the trace for anything in this plugin).
            if trace < 1e-10:
                self.qsi.debug('renormalizing state density matrix, trace was {}...', trace)
                self.qsi.sdm.renormalize()

            # The qubit is now no longer relevant in the SDM, at least until
            # the next gate is applied. So we can take it out.
            self.qs_ref = None

        # If this qubit was "classical", still process the method.
        else:
            if method == 'random':
                pass
            elif method == 'probable':
                pass
            elif method == 0:
                if self.classical != 0:
                    raise ValueError('cannot project to 0, qubit state was classical 1')
            elif method == 1:
                if self.classical != 1:
                    raise ValueError('cannot project to 1, qubit state was classical 0')
            else:
                raise ValueError('unknown measurement method %r' % method)

            # The probability is always 1.0
            p = 1.0

        return Measurement(self.qubit_ref, self.classical, struct.pack('<d', p), probability=p)

    def prep(self):
        """Put this qubit in the |0> state."""

        # Reset idle time.
        self._idle_time = 0

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

    def apply_pending_error(self):
        """Apply any remaining idling errors."""

        # If no time has passed, we don't have to do anything.
        if not self._idle_time:
            return

        # If this qubit is perfect, we don't have to do anything.
        if self.t1 == self.qsi.np.inf and self.t_phi == self.qsi.np.inf:
            return

        # Make sure the qubit is in the SDM.
        self.ensure_in_sdm()

        # Construct the Pauli transfer matrix for the error.
        gamma = 1 - self.qsi.np.exp(-self._idle_time / self.t1)
        lamda = 1 - self.qsi.np.exp(-self._idle_time / self.t_phi)
        ptm = self.qsi.ptm.amp_ph_damping_ptm(gamma, lamda)

        # Apply the Pauli transfer matrix.
        self.qsi.sdm.apply_ptm(self.qs_ref, ptm)

        # Reset the idle time.
        self._idle_time = 0

    def idle(self, time):
        """Pend some idle time for this qubit."""
        self._idle_time += time
