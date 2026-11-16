# SPDX-License-Identifier: Apache-2.0
"""python -m qoptct.estimate --toffoli T --qubits Q [--depth D] [--model ...] [--json] | --check"""
import sys

from .cli import main_estimate

sys.exit(main_estimate())
