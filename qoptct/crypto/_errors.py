# SPDX-License-Identifier: Apache-2.0
"""Exceptions shared by the qoptct.crypto generators."""


class CheckError(Exception):
    """A violated AND / erase precondition or a wrong value found by the classical (bit-sliced) run of a generated
    circuit (KaliskiCircuit.evaluate, MulSqCircuit.simulate)."""
