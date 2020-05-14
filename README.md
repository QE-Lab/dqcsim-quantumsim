# DQCsim backend wrapper for QuantumSim

[![PyPi](https://badgen.net/pypi/v/dqcsim-quantumsim)](https://pypi.org/project/dqcsim-quantumsim/)

See also: [DQCsim](https://github.com/mbrobbel/dqcsim) and
[QuantumSim](https://gitlab.com/quantumsim/quantumsim).

This repository contains some glue code to connect DQCsim to QuantumSim 0.2.0.

## Status

Roughly the following things are supported, but not tested for correctness or
reviewed by someone who knows what they're doing quantum-mechanics-wise:

 - Qubit amplitude damping and phase damping can be modelled and configured
   using t1/t2 times in DQCsim cycles on a per-qubit basis.
 - Measurements support collapsing the state randomly (default), selecting the
   most probable outcome, or collapsing to a predetermined state. The
   probability is reported either way.
 - All one- and two-qubit gates should be supported.

Not supported/limitations:

 - Modelling of gate-based errors.
 - Multi-qubit gates with more than two qubits (the QuantumSim API calls for
   this are missing).
 - Measurement gates affecting multiple qubits still collapse the state one at
   a time; it is currently not possible to get the full probability matrix for
   collapsing multiple states at once, and the `probable` measurement method
   also operates on a per-qubit basis.

## Install

You can install using `pip` using `pip install dqcsim-quantumsim` or
equivalent. If you're installing with `--user`, make sure that the path Python
installs the executables into is in your system path, otherwise DQCsim will not
be able to find the plugin. A simple way to see where the files are installed
is to run `pip uninstall dqcsim-quantumsim`; it shows which files it's about to
delete before querying for confirmation.

## Building/installing from source

 - "Build" the wheel file locally (necessary because of the executable shell
   script needed for DQCsim to recognize the plugin):
   `python3 setup.py bdist_wheel`

 - Install the wheel file you just built:
   `pip install target/python/dist/*` (or equivalent)

It should be safe to release the generated wheel to PyPI. Before doing so, make
sure to increment the version number in both `setup.py` and
`dqcsim_quantumsim/backend.py`.

## Usage

Once the plugin is installed, you can use it by selecting the `quantumsim`
backend in a DQCsim command.

The t1/t2 times of the qubits can be configured using an init arb (setting the
default for all qubits) or using an arb when allocating the qubit. The
interface/operation pair is `quantumsim.error`. The arb should have a JSON
object attached to it of the form `{"t1": <float>, "t2": <float>}`.

To control the way superposition is collapsed when a qubit is measured, attach
a JSON object of the form `{"method": ...}` to it, where method may be one of
the following:

 - `"random"`: collapse the state randomly based on the probabilities in the
   state vector.
 - `"probable"`: collapse to the most probable outcome.
 - an integer: the integer is decoded to binary representation, with the first
   qubit in the measurement operation mapping to the least significant bit. The
   state is always collapsed to this outcome.
 - a list: controls the method on a per-qubit basis when multiple qubits are
   measured at once. The list must be the same length as the number of qubits.
   It may be one of the string methods, the integer 0, or the integer 1.

If the probability for collapsing to a predetermined state is 0, the simulation
crashes. The probability for the selected measurement outcome is stored along
with the measurement result through a JSON object in the ArbData of the form
`{"probability": <float>}` and through the first binary string of the ArbData
as a C double.

## License

License is GPL, since QuantumSim's license is GPL.
