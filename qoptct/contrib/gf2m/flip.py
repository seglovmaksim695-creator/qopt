# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""Flip graph search (Kauers-Moosbauer) for low-rank decompositions of a bilinear map over GF(2).

A term is (a, b, c) int masks; the tensor is sum_r a_r (x) b_r (x) c_r.  A flip rewrites two terms sharing a factor
without changing the sum or the rank; a reduction merges two terms sharing two factors (rank - 1); a plus
transition splits a term (rank + 1) to leave a plateau.  Pure Python, intended for small tensors.
"""
import random
import time


def _reduce(D):
    """merge terms sharing two factors, drop zero terms; returns a new list."""
    changed = True
    D = [t for t in D if t[0] and t[1] and t[2]]
    while changed:
        changed = False
        seen = {}
        for idx, (a, b, c) in enumerate(D):
            for key, pos in (((0, a, b), 2), ((1, a, c), 1), ((2, b, c), 0)):
                if key in seen:
                    j = seen[key]
                    t1, t2 = list(D[j]), D[idx]
                    t1[pos] ^= t2[pos]
                    D[j] = tuple(t1)
                    D.pop(idx)
                    D = [t for t in D if t[0] and t[1] and t[2]]
                    changed = True
                    break
                seen[key] = idx
            if changed:
                break
    return D


def weight(D):
    return sum(bin(x).count("1") for t in D for x in t)


def search(dec, deadline, seed=1, plateau=20000, target=0):
    """random walk on the flip graph until `deadline` (time.time()); returns the best decomposition found
    (minimum rank, then minimum weight)."""
    rng = random.Random(seed)
    D = _reduce(list(dec))
    best = list(D)
    since = 0
    steps = 0
    while time.time() < deadline and len(best) > target:
        steps += 1
        R = len(D)
        if R < 2:
            break
        k = rng.randrange(3)
        groups = {}
        for i, t in enumerate(D):
            groups.setdefault(t[k], []).append(i)
        pairs = [g for g in groups.values() if len(g) > 1]
        if not pairs or since > plateau:
            # plus transition
            i, j = rng.sample(range(R), 2)
            a1, b1, c1 = D[i]
            a2 = D[j][0]
            if a1 != a2:
                D[i] = (a1 ^ a2, b1, c1)
                D.append((a2, b1, c1))
            since = 0
            D = _reduce(D)
            continue
        g = rng.choice(pairs)
        i, j = rng.sample(g, 2)
        t1, t2 = list(D[i]), list(D[j])
        o1, o2 = [x for x in range(3) if x != k]
        if rng.random() < 0.5:
            o1, o2 = o2, o1
        # t1 = (s, u1, v1), t2 = (s, u2, v2) -> (s, u1 + u2, v1), (s, u2, v2 + v1)
        t1[o1] ^= t2[o1]
        t2[o2] ^= t1[o2]
        D[i], D[j] = tuple(t1), tuple(t2)
        if not (t1[o1] and t2[o2]):
            D = _reduce(D)
        elif steps % 8 == 0:
            D = _reduce(D)
        since += 1
        if len(D) < len(best) or (len(D) == len(best) and weight(D) < weight(best)):
            if len(D) < len(best):
                since = 0
            best = list(D)
    return best
