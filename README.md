# DQCsim backend wrapper for QuantumSim

[![PyPi](https://badgen.net/pypi/v/dqcsim-quantumsim)](https://pypi.org/project/dqcsim-quantumsim/)

See also: [DQCsim](https://github.com/mbrobbel/dqcsim) and
[QuantumSim](https://gitlab.com/quantumsim/quantumsim).

This is VERY alpha right now, built for DQCsim 0.0.2 and QuantumSim 0.2.0 and
hardly tested, though preliminary testing shows that it seems to work.

No error modelling is supported yet. This is just a dumb wrapper that can
execute one- and two-qubit gates received from the upstream plugin using
QuantumSim.

USE AT YOUR OWN RISK.

License is GPL, since QuantumSim's license is GPL.
