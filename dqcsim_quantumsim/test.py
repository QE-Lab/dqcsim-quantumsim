import unittest
from dqcsim.plugin import *
from dqcsim.host import *
import tempfile
import os

@plugin("Deutsch-Jozsa", "Tutorial", "0.1")
class DeutschJozsa(Frontend):

    def oracle_constant_0(self, qi, qo):
        """x -> 0 oracle function."""
        pass

    def oracle_constant_1(self, qi, qo):
        """x -> 1 oracle function."""
        self.x_gate(qo)

    def oracle_passthrough(self, qi, qo):
        """x -> x oracle function."""
        self.cnot_gate(qi, qo)

    def oracle_invert(self, qi, qo):
        """x -> !x oracle function."""
        self.cnot_gate(qi, qo)
        self.x_gate(qo)

    def deutsch_jozsa(self, qi, qo, oracle, expected):
        """Runs the Deutsch-Jozsa algorithm on the given oracle. The oracle is
        called with the input and output qubits as positional arguments."""

        # Prepare the input qubit.
        self.prepare(qi)
        self.h_gate(qi)

        # Prepare the output qubit.
        self.prepare(qo)
        self.x_gate(qo)
        self.h_gate(qo)

        # Run the oracle function.
        oracle(qi, qo)

        # Measure the input.
        self.h_gate(qi)
        self.measure(qi)
        measurement = self.get_measurement(qi)
        if measurement.value:
            self.info('Oracle was balanced! The probability is %f' % measurement['probability'])
            if expected != 'balanced':
                raise ValueError('unexpected oracle result!')
        else:
            self.info('Oracle was constant! The probability is %f' % measurement['probability'])
            if expected != 'constant':
                raise ValueError('unexpected oracle result!')

    def handle_run(self):
        qi, qo = self.allocate(2)

        self.info('Running Deutsch-Jozsa on x -> 0...')
        self.deutsch_jozsa(qi, qo, self.oracle_constant_0, 'constant')

        self.info('Running Deutsch-Jozsa on x -> 1...')
        self.deutsch_jozsa(qi, qo, self.oracle_constant_1, 'constant')

        self.info('Running Deutsch-Jozsa on x -> x...')
        self.deutsch_jozsa(qi, qo, self.oracle_passthrough, 'balanced')

        self.info('Running Deutsch-Jozsa on x -> !x...')
        self.deutsch_jozsa(qi, qo, self.oracle_invert, 'balanced')

        self.free(qi, qo)


class Constructor(unittest.TestCase):

    def test_simple(self):
        with Simulator(
            (DeutschJozsa(), {
                'verbosity': Loglevel.INFO
            }),
            ('quantumsim', {
                'verbosity': Loglevel.INFO
            }),
            stderr_verbosity=Loglevel.INFO
        ) as sim:
            sim.run()
