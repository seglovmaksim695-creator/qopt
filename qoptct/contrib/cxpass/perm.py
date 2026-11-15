# SPDX-License-Identifier: Apache-2.0
"""Absorption of an ancilla permutation: A' = P A, the permutation is carried over to all later gates."""
import random


def best_perm(A, nd, rng=None, tries=1):
    """A: rows (masks over local wire indices). Returns pi: dict old wire -> new wire (only for wires >= nd) such that
    the rows A'[pi[w]] = A[w] have as many ones on the diagonal as possible."""
    n = len(A)
    anc = [w for w in range(n) if w >= nd]
    rng = rng or random.Random(0)
    adj = {w: [x for x in range(nd, n) if A[w] >> x & 1] for w in anc}
    match = {}  # x -> w
    def aug(w, seen):
        order = adj[w][:]; rng.shuffle(order)
        if w in order: order.remove(w); order.insert(0, w)
        for x in order:
            if x in seen: continue
            seen.add(x)
            if x not in match or aug(match[x], seen):
                match[x] = w; return True
        return False
    for w in sorted(anc, key=lambda w: len(adj[w])):
        aug(w, set())
    pi = {}
    for x, w in match.items(): pi[w] = x
    free = [x for x in anc if x not in match]
    rest = [w for w in anc if w not in pi]
    for w in rest[:]:
        if w in free: pi[w] = w; free.remove(w); rest.remove(w)
    for w, x in zip(rest, free): pi[w] = x
    return pi


def apply_perm_rows(A, pi, nd):
    n = len(A)
    A2 = list(A)
    for w in range(n):
        x = pi.get(w, w)
        A2[x] = A[w]
    return A2
