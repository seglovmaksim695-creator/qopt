"""check_bloq: prove a Qualtran bloq against Qualtran's own classical action.

What "PROVED" means here
------------------------
The bloq is lowered gate by gate with the strict bridge (rules A1 controlled-And and C1
QROAMCleanAdjoint are enabled), translated to the logical language of qoptct
(``x``, ``cx``, ``and``, ``erase``) and then two independent checks are made on the WHOLE
classical input domain of the bloq (all values of all LEFT/THRU registers, or the explicit
``domain`` you pass; the classical domain of a register dtype, e.g. ``1..mod-1`` for a Montgomery type, also
restricts the claim and is reported in the result field ``domain``):

1. exhaustive simulation of the program on every input, with every measurement outcome
   (an ``erase`` with a non-zero phase coefficient would give a sign on some outcome, so ANY
   non-zero coefficient fails) -- the final state must be exactly the output of
   ``bloq.call_classically`` on the output wires, and ZERO on every other wire;
2. the qoptct BDD prover ``kprove`` (unchanged) on the same program with a specification
   built from the same classical table, plus the exact partial-domain guard of the bridge
   (every And target is clean and every erase condition holds exactly).

Verdicts: ``PROVED``, ``REFUTED`` (with a witness input) and ``UNKNOWN`` (the bridge does not
support some part, ``call_classically`` is not defined on part of the domain, the domain is
too large, or a resource budget was hit).  ``UNKNOWN`` is never a proof.

Assumptions that are recorded in the result (``assumptions``)
-------------------------------------------------------------
* Qualtran's ``call_classically`` is the specification (we do not audit it);
* the rules of the bridge (``bridge/rules.py``, A1 ``controlled_and``, C1 ``qroam_clean_adjoint``);
* the unitary action on superpositions follows from the action on basis states because the
  program is built only from classical reversible gates and X-basis measurements with explicit
  phase corrections (the simulator reports the sign of every outcome branch).

Only small sizes are meant: the number of input bits is limited by ``max_input_bits``.
"""
from __future__ import annotations

import itertools
from hashlib import sha256

__all__ = ['check_bloq']


def _versions():
    out = {}
    try:
        from importlib.metadata import version
        for name in ('qualtran', 'z3-solver', 'numpy'):
            try:
                out[name] = version(name)
            except Exception:
                out[name] = None
    except Exception:  # pragma: no cover
        pass
    return out


def _bits_msb(dtype, value):
    bits = list(dtype.to_bits(value))
    return [int(b) for b in bits]


def check_bloq(bloq, *, domain=None, max_input_bits=12, use_prover=True, budget=120, max_mb=900,
               extra_atom_rule=None, max_leaves=200000):
    """Prove ``bloq`` against ``bloq.call_classically`` on its whole input domain.

    Parameters
    ----------
    bloq : qualtran.Bloq
        A bloq with only pure quantum registers (no classical bits).  Shaped registers are allowed.
    domain : callable or None
        Optional promise ``domain(values) -> bool`` where ``values`` maps each LEFT/THRU register
        name to its classical value (an int, or a numpy object array for shaped registers).
        Inputs for which it returns False are excluded from the claim; the claim is then
        "for all inputs satisfying the promise" and the result says so.
    max_input_bits : int
        Refuse (UNKNOWN) when the total number of input bits is larger.
    use_prover : bool
        Also run the qoptct BDD prover (recommended; a PROVED verdict needs both checks).
        When False the verdict is capped at ``UNKNOWN`` for a passing enumeration.
    budget, max_mb : prover time (s) and memory (MB) limits.
    extra_atom_rule : callable or None
        Atom rule for the bridge; default is rule A1 + rule C1.

    Returns
    -------
    dict with ``verdict`` in {'PROVED', 'REFUTED', 'UNKNOWN'}, ``detail``, ``witness`` (REFUTED),
    ``inputs_checked``, ``qasm_sha256``, ``layout_sha256``, ``prover`` (kprove result),
    ``assumptions`` and ``versions``.
    """
    import numpy as np
    from qualtran import Side

    from .bridge.from_qualtran import lower_graph
    from .bridge.qroam_clean_adjoint import chained_atom_rule
    from .bridge.rules import TranslationError, translate
    from .bridge.sim import SimulationError, run
    from .bridge.strict_proof import digest_json, prove_program

    result = dict(verdict='UNKNOWN', detail='', witness=None, inputs_checked=0, method=[],
                  assumptions=['Qualtran call_classically is the specification',
                               'bridge rules (A1 controlled-And, C1 QROAMCleanAdjoint) are trusted as documented',
                               'action on superpositions follows from basis states plus the phase check of all outcomes'],
                  versions=_versions(), bloq=repr(bloq)[:200],
                  domain='all bit patterns of the inputs' if domain is None else 'inputs satisfying the given promise')

    def unknown(detail):
        result.update(verdict='UNKNOWN', detail=detail)
        return result

    # --- lowering and translation -------------------------------------------------------------------
    try:
        lowered = lower_graph(bloq, max_leaves=max_leaves, extra_atom_rule=extra_atom_rule or chained_atom_rule())
        prog = translate(lowered, lowered['inputs'], lowered['outputs'])
    except TranslationError as exc:  # LayoutError is a TranslationError
        return unknown(f'bridge refused the bloq: {type(exc).__name__}: {exc}')
    except Exception as exc:
        return unknown(f'lowering failed: {type(exc).__name__}: {exc}')
    result.update(nq=prog['nq'], ops=len(prog['ops']), qasm_sha256=sha256(prog['qasm'].encode()).hexdigest(),
                  layout_sha256=digest_json(prog['layout']))

    in_regs = [r for r in lowered['registers']['inputs']]
    out_regs = [r for r in lowered['registers']['outputs']]
    lefts = {r.name: r for r in bloq.signature.lefts()}
    rights = {r.name: r for r in bloq.signature.rights()}
    total_bits = sum(r['bitsize'] * len(r['elements']) for r in in_regs)
    result['input_bits'] = total_bits
    if total_bits > max_input_bits:
        return unknown(f'{total_bits} input bits > max_input_bits={max_input_bits}')
    for reg in list(lefts.values()) + list(rights.values()):
        if reg.dtype.num_cbits != 0:
            return unknown('classical bits in the signature are not supported by check_bloq')

    # --- the classical domain -----------------------------------------------------------------------
    elements = []  # (register name, index tuple, dtype, wires msb)
    for desc in in_regs:
        reg = lefts[desc['name']]
        for el in desc['elements']:
            elements.append((desc['name'], tuple(el['index']), reg.dtype, el['wires_msb']))
    try:
        element_domains = [list(dt.get_classical_domain()) for _, _, dt, _ in elements]
    except Exception as exc:
        return unknown(f'no classical domain for an input register: {type(exc).__name__}: {exc}')
    count = 1
    for d in element_domains:
        count *= len(d)
    if count > (1 << max_input_bits):
        return unknown(f'{count} classical inputs > 2^{max_input_bits}')

    out_elements = []
    for desc in out_regs:
        reg = rights[desc['name']]
        for el in desc['elements']:
            out_elements.append((desc['name'], tuple(el['index']), reg.dtype, el['wires_msb']))

    def to_values(assignment):
        vals = {}
        for (name, idx, _dt, _w), value in zip(elements, assignment):
            if lefts[name].shape:
                if name not in vals:
                    vals[name] = np.empty(lefts[name].shape, dtype=object)
                vals[name][idx] = value
            else:
                vals[name] = value
        return vals

    nq = prog['nq']
    in_wire_order = []  # one entry per input wire: wire index; position defines the table index bit
    for _, _, dt, wires in elements:
        in_wire_order.extend(reversed(wires))  # LSB first within each element
    # table[wire] = list over assignments (in enumeration order) of the expected output bit
    expected_bits_all = []
    inputs_checked = 0
    rows = []  # (initial state int, expected final int, assignment dict for the witness)
    total_assignments = 0
    for assignment in itertools.product(*element_domains):
        total_assignments += 1
        vals = to_values(assignment)
        if domain is not None:
            try:
                if not domain(vals):
                    continue
            except Exception as exc:
                return unknown(f'domain promise raised {type(exc).__name__}: {exc}')
        try:
            classical = bloq.call_classically(**vals)
            if not isinstance(classical, tuple):
                classical = (classical,)
            names = [r.name for r in bloq.signature.rights()]
            if len(classical) != len(names):
                return unknown(f'call_classically returned {len(classical)} values for {len(names)} output registers')
            outs = dict(zip(names, classical))
        except Exception as exc:
            return unknown('call_classically is not defined on the whole domain (use `domain=` to state a promise); '
                           f'first failure at {_show(vals)}: {type(exc).__name__}: {exc}')
        initial = 0
        for (name, idx, dt, wires), value in zip(elements, assignment):
            bits = _bits_msb(dt, value)
            if len(bits) != len(wires):
                return unknown(f'register {name}: unexpected bit width of a classical value')
            for bit, wire in zip(bits, wires):
                initial |= bit << wire
        final = 0
        try:
            for name, idx, dt, wires in out_elements:
                value = outs[name][idx] if rights[name].shape else outs[name]
                bits = _bits_msb(dt, value)
                if len(bits) != len(wires):
                    return unknown(f'output {name}: unexpected bit width of a classical value')
                for bit, wire in zip(bits, wires):
                    final |= bit << wire
        except Exception as exc:
            return unknown(f'call_classically output is not encodable at {_show(vals)}: {type(exc).__name__}: {exc}')
        rows.append((initial, final, vals))
    if not rows:
        return unknown('the promise excludes every input')

    # --- check 1: exhaustive simulation with phases of every outcome ---------------------------------
    for initial, final, vals in rows:
        try:
            got, phases = run(nq, prog['ops'], initial)
        except SimulationError as exc:
            result.update(verdict='REFUTED', witness=_show(vals), detail=f'program is not executable on this input: {exc}',
                          inputs_checked=inputs_checked)
            return result
        if got != final or any(phases):
            result.update(verdict='REFUTED', witness=_show(vals), inputs_checked=inputs_checked,
                          detail=('wrong final state (expected/got as wire masks): '
                                  f'{final:#x} / {got:#x}' if got != final else
                                  f'uncorrected measurement phase on outcome branches: {phases}'))
            return result
        inputs_checked += 1
    result['inputs_checked'] = inputs_checked
    if len(rows) != (1 << total_bits):
        result['domain'] = (f'{len(rows)} of {1 << total_bits} bit patterns: the classical domain of the register '
                            'dtypes' + ('' if domain is None else ' restricted by the given promise'))
    result['method'].append('exhaustive_simulation_all_outcomes_with_phase')

    if not use_prover:
        return unknown('enumeration passed; prover not requested (verdict stays UNKNOWN)')

    # --- check 2: the qoptct BDD prover on the same program -----------------------------------------
    from qoptct.verify.bdd import ONE, ZERO

    initial_to_final = {}
    for initial, final, _ in rows:
        initial_to_final[initial] = final
    layout = prog['layout']
    in_wires = list(in_wire_order)
    out_wires = sorted({q for word in layout['outputs'].values() for q in word})

    def bit_table(wire):
        # index bit j of an assignment corresponds to in_wires[j]
        table = [0] * (1 << len(in_wires))
        for initial, final in initial_to_final.items():
            index = 0
            for j, w in enumerate(in_wires):
                index |= ((initial >> w) & 1) << j
            table[index] = (final >> wire) & 1
        return table

    def build(S, table, k):
        # Shannon expansion on the highest remaining wire
        if k == 0:
            return ONE if table[0] else ZERO
        half = 1 << (k - 1)
        hi = build(S, table[half:], k - 1)
        lo = build(S, table[:half], k - 1)
        if hi == lo:
            return hi
        return S.B.ITE(S.x(in_wires[k - 1]), hi, lo)

    def symbolic(S):
        return {q: build(S, bit_table(q), len(in_wires)) for q in out_wires}

    pre = None
    if len(rows) != (1 << len(in_wires)):
        allowed = [0] * (1 << len(in_wires))
        for initial in initial_to_final:
            index = 0
            for j, w in enumerate(in_wires):
                index |= ((initial >> w) & 1) << j
            allowed[index] = 1

        def pre(S):
            return build(S, allowed, len(in_wires))

    contract = dict(source='qualtran call_classically', bloq=repr(bloq)[:200],
                    domain=result['domain'], input_bits=total_bits, inputs=len(rows),
                    table_sha256=digest_json(sorted(initial_to_final.items())))
    receipt = dict(status='PASS', qasm_sha256=result['qasm_sha256'], contract_sha256=digest_json(contract),
                   layout_sha256=result['layout_sha256'], independent_inputs=inputs_checked,
                   all_wires_checked=True, phase_checked=True)
    proof = prove_program(prog['qasm'], layout, symbolic, pre=pre, validation=receipt, contract=contract,
                          budget=budget, max_mb=max_mb)
    result['prover'] = {k: v for k, v in proof.items() if k in
                        ('status', 'detail', 'prover_sha256', 'wall_seconds', 'peak_rss_mb', 'method', 'artifact',
                         'prover_invoked', 'witness')}
    result['method'].append('qoptct_kprove_bdd_with_exact_partial_domains')
    if proof.get('partial_domains'):
        result['partial_domains'] = {k: proof['partial_domains'].get(k) for k in
                                     ('status', 'domain_checks', 'bdd_nodes')}
    status = proof.get('status')
    if status == 'PROVED':
        result.update(verdict='PROVED', detail='exhaustive simulation and the BDD prover agree')
    elif status == 'REFUTED':
        result.update(verdict='REFUTED', witness=proof.get('witness'), detail='the BDD prover refuted the program')
    else:
        result.update(verdict='UNKNOWN', detail=f'prover status {status}: {proof.get("detail")}')
    return result


def _show(vals):
    out = {}
    for k, v in vals.items():
        out[k] = v.tolist() if hasattr(v, 'tolist') else v
    return out
