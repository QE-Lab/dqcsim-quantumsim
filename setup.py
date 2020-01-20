#!/usr/bin/env python3

import os
from setuptools import setup

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()

setup(
    name = "dqcsim-quantumsim",
    version = "0.0.3",
    author = "Jeroen van Straten",
    author_email = "j.vanstraten-1@tudelft.nl",
    description = "DQCsim backend for QuantumSim.",
    license = "GPLv3",
    keywords = "dqcsim quantumsim",
    url = "https://github.com/QE-Lab/dqcsim-quantumsim",
    long_description = read('README.md'),
    long_description_content_type = 'text/markdown',
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering",
    ],
    packages = ['dqcsim_quantumsim'],
    install_requires = ['dqcsim>=0.0.13', 'quantumsim==0.2'],
    tests_require = [
        'nose',
    ],
    test_suite = 'nose.collector',
    data_files = [
        ('bin', [
            'data/bin/dqcsbequantumsim',
        ]),
    ],
)
