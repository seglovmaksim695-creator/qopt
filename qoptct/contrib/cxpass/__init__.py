# SPDX-License-Identifier: Apache-2.0
"""cxpass -- CX-reduction post-passes for OpenQASM 3 circuits with measurements and conditional corrections.

    improve(qasm, clean, budget_sec=60.0, ref=None, n_data=None) -> str

Input: the OpenQASM 3 subset written by qoptct (``--mode t`` / ``balanced`` / the c2, arith, gf2m compilers): gates
h x z s sdg t tdg cx cz, ``c[k] = measure q[j];``, ``reset q[j];`` and ``if (cond) { ... }`` with AND/XOR conditions.
Output: a text that is never worse in (T, CX) -- both counts are <= the input and at least one is smaller -- and that passed
verification; otherwise the input text is returned unchanged (also on any error or when the budget runs out).

Passes (from the qoptct lab, see lab/v04/cxpass.md):
  * resynth -- block resynthesis of CX/CZ/phase regions with recomputed Clifford weights, ancilla-permutation absorption and
               |+>-wire freedom (agent b10);
  * merge   -- conditional CZ on the same pair of wire forms merged into one CZ with an XOR condition (agent a7);
  * cpo     -- conditional phase-correction optimiser: corrections carried forward as polynomials of path variables and emitted
               together with XOR conditions (agent a10);
  * etr     -- E-transform of the final CX block before X-measured garbage (agent a7);
  * corr    -- conditional CZ = S S Sdg on wire parities (agent a7).

Verification: with ``ref`` (path of the .qc reference) the result is checked against it in the strict (clean=False) or clean
model (qoptct.verify.prove when installed, else the state-vector checker of cxpass.check, else for strict circuits with more
than 22 qubits the c10 path-sum verifier).  Without ``ref`` the result is checked against the input program itself on the wires
on which the input is deterministic.  Results may contain XOR conditions (``if (c[0] ^ c[3]) {...}``); verifiers that read
them: cxpass.check, qoptct.contrib.c10.verify_big, qoptct.verify.
"""
import time

from . import ir

__all__ = ["improve", "improve_report", "PASSES"]

PASSES = ("resynth", "merge", "cpo", "etr", "corr")


def _key(text):
    return ir.counts_text(text)


def _not_worse(new, old):
    return new[0] <= old[0] and new[1] <= old[1]


def _better(new, old):
    return _not_worse(new, old) and new != old


def _run_pass(name, text, n_data, deadline, seed=1):
    """One pass on normalised text -> new text (not verified)."""
    N, items = ir.parse(text)
    if name == "merge":
        from .merge import merge
        return ir.emit(N, merge(items, N))
    if name == "corr":
        from .corr import relocate
        out, _ = relocate(list(items), N)
        return ir.emit(N, out)
    if name == "etr":
        from . import etr
        left = max(1.0, deadline - time.time())
        out = etr.optimize(list(items), N, seed=seed, rounds=4, tlimit=left)
        return ir.emit(N, out)
    if name == "cpo":
        from . import cpo
        txt, _ = cpo.optimize(text, n=n_data)
        return ir.normalize(txt)
    if name == "resynth":
        from . import resynth, blocks
        head, nq, out, _ = resynth.process(text, n_data if n_data is not None else N, seed=seed, tries=6, deadline=deadline)
        return ir.normalize(blocks.emit_qasm(head, nq, out))
    raise ValueError(name)


def improve_report(qasm, clean=False, budget_sec=60.0, ref=None, n_data=None, passes=PASSES, rounds=3, log=None,
                   official_verify=True, allow_xor=True):
    """Like improve(), but returns (text, report dict)."""
    t0 = time.time(); total = max(1.0, float(budget_sec))
    deadline = t0 + 0.7 * total          # search; the rest is kept for the final check
    hard = t0 + 0.95 * total
    rep = {"in": None, "out": None, "steps": [], "verified": False}
    try:
        k_in = _key(qasm); rep["in"] = k_in; rep["out"] = k_in
        N, items = ir.parse(qasm)
        cur = ir.emit(N, items)
        if ref is not None:
            import qoptct
            n_data = len(qoptct.load_qc(ref)[0])
        D = None
        if ref is None:
            from . import check
            if N > check.MAX_SV:
                rep["error"] = "no reference and too many qubits to check against the input"
                return qasm, rep
            cand = list(range(n_data)) if (n_data is not None and not clean) else list(range(N))
            D = check.deterministic_wires(N, items, cand, greedy=clean, deadline=t0 + 0.3 * (deadline - t0))
            if not D:
                rep["error"] = "input program is not deterministic"
                return qasm, rep
            if n_data is not None and not clean and D != list(range(n_data)):
                rep["error"] = "input program is not deterministic on its data wires"
                return qasm, rep
        from . import check

        def ok(text, budget):
            return check.verify(text, clean, ref_qc=ref, ref_text=None if ref else cur0, D=D, budget=budget,
                                official=official_verify)
        cur0 = cur
        k_cur = _key(cur)
        dead = set()   # passes that produced nothing (or failed) on the current text
        for r in range(rounds * len(passes)):
            # every pass is run on the current text; the best candidate that verifies is taken (greedy, best first)
            cands = []
            for p in passes:
                if p in dead or time.time() > deadline:
                    continue
                t1 = time.time()
                try:
                    new = _run_pass(p, cur, n_data, deadline, seed=1 + r)
                    if ir.parse(new)[0] != N:
                        raise ValueError("qubit count changed")
                    if not allow_xor and "^" in new and "^" not in qasm:
                        raise ValueError("XOR condition not allowed")
                except Exception as e:  # a pass that fails is skipped
                    rep["steps"].append({"pass": p, "error": repr(e)[:200]}); dead.add(p)
                    continue
                k_new = _key(new)
                step = {"pass": p, "before": k_cur, "after": k_new, "sec": round(time.time() - t1, 2)}
                rep["steps"].append(step)
                if log:
                    log(step)
                if _better(k_new, k_cur):
                    cands.append((k_new, p, new, step))
                else:
                    dead.add(p)
            improved = False
            for k_new, p, new, step in sorted(cands, key=lambda c: (c[0][1], c[0][0])):
                if time.time() > deadline:
                    break
                good = ok(new, max(2.0, min(120.0, hard - time.time())))
                step["verified"] = bool(good)
                if good:
                    cur, k_cur, improved = new, k_new, True
                    dead = set()
                    break
                dead.add(p)
            if not improved:
                break
        if _better(k_cur, k_in):
            # final, stronger check of the chosen text (more random states / measurement branches)
            if not check.verify(cur, clean, ref_qc=ref, ref_text=None if ref else cur0, D=D,
                                budget=max(1.0, min(120.0, hard - time.time())),
                                official=official_verify, trials=30, strong=True):
                rep["error"] = "final verification failed"
                rep["sec"] = round(time.time() - t0, 2)
                return qasm, rep
            rep["out"] = k_cur; rep["verified"] = True
            rep["sec"] = round(time.time() - t0, 2)
            return cur, rep
    except Exception as e:
        rep["error"] = repr(e)[:300]
    rep["sec"] = round(time.time() - t0, 2)
    return qasm, rep


def improve(qasm: str, clean: bool, budget_sec: float = 60.0, ref=None, n_data=None, allow_xor=True) -> str:
    """CX post-pass. Returns a verified text with T and CX not above the input (at least one smaller), else ``qasm`` itself.

    ref:    optional path of the .qc reference (recommended: verification in the requested model; enables ancilla tricks).
    n_data: optional number of data wires (qubits 0..n_data-1); taken from ``ref`` when it is given.
    allow_xor: False forbids XOR conditions in the result (for consumers whose checker only reads AND conditions).
    """
    try:
        return improve_report(qasm, clean, budget_sec, ref=ref, n_data=n_data, allow_xor=allow_xor)[0]
    except Exception:
        return qasm
