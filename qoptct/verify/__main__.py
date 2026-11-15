# SPDX-License-Identifier: Apache-2.0
"""python -m qoptct.verify CIRCUIT.qasm REFERENCE.qc [--clean] [--budget SEC] [--json]"""
import sys

from .cli import main_prove

sys.exit(main_prove())
