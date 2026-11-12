"""Sparse application of the project's one- and two-qubit gates.

Each matrix row has at most two nonzero entries for the current gate set.
This avoids multiplying through the zeros of the dense matrices.
"""

import numpy as np


class FastGateBatch:
    def __init__(self, matrices):
        matrices = np.asarray(matrices)
        if matrices.ndim != 3 or matrices.shape[1] != matrices.shape[2]:
            raise ValueError("Expected a batch of square gate matrices")
        counts = np.count_nonzero(matrices, axis=2)
        if np.any(counts > 2):
            raise ValueError("FastGateBatch supports at most two nonzeros per row")
        self.dim = matrices.shape[1]
        self.indices = np.zeros((*matrices.shape[:2], 2), dtype=np.intp)
        self.values = np.zeros((*matrices.shape[:2], 2), dtype=matrices.dtype)
        for a, matrix in enumerate(matrices):
            for row in range(self.dim):
                columns = np.flatnonzero(matrix[row])
                self.indices[a, row, :len(columns)] = columns
                self.values[a, row, :len(columns)] = matrix[row, columns]

    def __call__(self, states):
        states = np.asarray(states)
        if states.ndim != 2 or states.shape[1] != self.dim:
            raise ValueError("Unexpected state batch shape")
        return np.einsum("ark,bark->bar", self.values, states[:, self.indices])
