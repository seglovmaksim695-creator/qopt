"""Targeted corruptions: exact rare inputs, erase corrections, segment joints.

Does not call the prover. A witness is obtained by executing the real mutant with
the independent simulator; an error of the simulator itself is not a witness.
The contract is passed explicitly, without assuming p=2**n-1. All output words,
explicit junk included, must be part of expected(values).
"""
from hashlib import sha256
from itertools import combinations

from .emit import emit_qasm
from .sim import run, SimulationError, parse_qasm


def encode(values, registers):
    result = 0
    for name, wires in registers.items():
        value = values[name]
        if type(value) is not int or not 0 <= value < 1 << len(wires):
            raise ValueError(f'Invalid word {name}: {value}')
        for bit, wire in enumerate(wires):
            result |= ((value >> bit) & 1) << wire
    return result


def decode(state, registers):
    return {name: sum(((state >> wire) & 1) << i for i, wire in enumerate(wires))
            for name, wires in registers.items()}


def generate(prog, expected, domain, samples, *, rare_cases=(), boundaries=(), count=20,
             correction_mode='remove',correction_groups=()):
    """rare_cases: [{label, match:{word:value}, witness:{all_inputs}}].

    boundaries: [{position, label, wires}]; position means just BEFORE op.
    At least half accepted mutants must target a boundary, correction or rare
    input. If this cannot be achieved, raise; never silently fill the quota.
    Each returned record contains actual raw QASM hashes and witness.
    """
    original = list(prog['ops'])
    layout = prog['layout']
    nq = prog['nq']
    baseline_text=prog['qasm']
    if baseline_text!=emit_qasm(nq,original):
        raise ValueError('The baseline QASM passed in does not match the exact IR')
    parsed_nq,baseline_ops=parse_qasm(baseline_text)
    if parsed_nq!=nq:
        raise ValueError('The independent parser found another size of the baseline')
    ins, outs = layout['inputs'], layout['outputs']
    if layout.get('nq')!=nq:
        raise ValueError('The layout size does not match')
    for registers in (ins,outs):
        flat=[q for word in registers.values() for q in word]
        if not flat or len(flat)!=len(set(flat)) or any(type(q) is not int or not 0<=q<nq for q in flat):
            raise ValueError('Overlapping or invalid wire of the layout')
    if layout.get('garbage'):
        raise ValueError('Arbitrary garbage is forbidden; give junk as a full output')
    if count < 20:
        raise ValueError('At least 20 canaries are required')
    if correction_mode not in ('remove','replace'):
        raise ValueError('Unknown way to corrupt corrections')
    output_wires = [q for word in outs.values() for q in word]
    allowed = sum(1 << q for q in output_wires)
    values_seen, prepared = set(), []
    for values in list(samples) + [r['witness'] for r in rare_cases]:
        if set(values) != set(ins):
            raise ValueError('Incomplete input of the independent check')
        key = tuple(values[k] for k in ins)
        if key in values_seen or not domain(values):
            continue
        values_seen.add(key)
        initial = encode(values, ins)
        want = expected(values)
        if set(want) != set(outs):
            raise ValueError('The contract does not define all outputs, junk included')
        final, phase = run(nq, baseline_ops, initial)
        if decode(final, outs) != want or final & ~allowed or any(phase):
            raise ValueError(f'The source circuit failed the independent check: {values}')
        prepared.append((dict(values), initial, want, final))
    if not prepared:
        raise ValueError('Empty input domain of the canaries')
    baseline = sha256(baseline_text.encode()).hexdigest()
    mutants, seen = [], {baseline}

    def accept(kind, changed, size, mutation, targeted, predicate=None):
        qasm = emit_qasm(size, changed)
        parsed_size,parsed_ops=parse_qasm(qasm)
        if parsed_size!=size:
            raise ValueError('The independent parser found another size of the mutant')
        digest = sha256(qasm.encode()).hexdigest()
        if digest in seen:
            return False
        seen.add(digest)
        bad, first_phase, phase_pair = None, None, None
        triggered, inactive = 0, 0
        try:
            for values, initial, want, basefinal in prepared:
                final, phase = run(size, parsed_ops, initial)
                got = decode(final, outs)
                if predicate is not None:
                    match, flag, target = predicate
                    trigger = all(values[k] == v for k,v in match.items())
                    delta = ((1 << flag) | (1 << target)) if trigger else 0
                    if final != basefinal ^ delta or any(phase):
                        return False
                    triggered += trigger
                    inactive += not trigger
                if bad is None and (got != want or final & ~allowed):
                    bad = dict(reason='wrong_output' if got != want else 'dirty_ancilla',
                               inputs=values,want=want,got=got,dirty_bits=final & ~allowed)
                ph = tuple(phase)
                phase_info = dict(inputs=values,measurements=len(phase),
                                  nonzero_coefficients=[i for i,v in enumerate(phase) if v])
                if first_phase is None:
                    first_phase = ph, phase_info
                elif ph != first_phase[0] and phase_pair is None:
                    phase_pair = dict(reason='input_dependent_phase',
                                      first=first_phase[1],second=phase_info)
        except SimulationError:
            return False
        if kind == 'erasure_correction':
            bad = phase_pair
        elif bad is None:
            bad = phase_pair
        if bad is None or (predicate is not None and (not triggered or not inactive)):
            return False
        witness = dict(bad,checked_inputs=len(prepared))
        if predicate is not None:
            witness.update(triggered_checked=triggered,untriggered_checked=inactive,
                           only_when=predicate[0])
        mutants.append(dict(kind=kind,targeted=targeted,nq=size,ops=changed,
                            witness=witness,mutation=mutation,qasm_sha256=digest,
                            baseline_qasm_sha256=baseline))
        return True

    # First preserve diversity: up to four different segment boundaries.
    for boundary in boundaries:
        pos = boundary['position']
        if type(pos) is not int or not 0 <= pos <= len(original):
            raise ValueError('Boundary outside the program')
        for q in boundary['wires']:
            if accept('segment_boundary',original[:pos]+[('x',q)]+original[pos:],nq,
                      dict(boundary=boundary['label'],position=pos,wire=q),True):
                break
        if len(mutants) >= 4:
            break

    # Explicit source groups each require their own independent phase witness.
    if correction_groups:
        from .correction_group_canaries import generate as grouped_corrections
        grouped_corrections(original,ins,nq,correction_groups,correction_mode,accept)
    # Default sampling preserves the original established control order.
    erasures = [] if correction_groups else [i for i,op in enumerate(original) if op[0]=='erase']
    if len(erasures)>40:
        erasures = [erasures[i*(len(erasures)-1)//39] for i in range(40)]
    phase_count = 0
    for i in erasures:
        _,t,corrections = original[i]
        for j in range(len(corrections)):
            if correction_mode=='remove':
                variants=[(tuple(corrections[:j])+tuple(corrections[j+1:]),
                           dict(position=i,removed=corrections[j]))]
            else:
                # Keep a canonical ERASE gadget for the unchanged Walker.
                # The wrong correction must be refuted semantically, not
                # merely refused because removal left its accepted grammar.
                variants=[]
                correction=corrections[j]
                occupied={t}|{q for item in corrections for q in item[1:]}
                for q in [q for word in ins.values() for q in word]:
                    if q in occupied:
                        continue
                    replacement=tuple(correction[:-1])+(q,)
                    variants.append((tuple(corrections[:j])+(replacement,)+tuple(corrections[j+1:]),
                                     dict(position=i,original=correction,replacement=replacement)))
                    if len(variants)>=12:
                        break
            for corrected,mutation in variants:
                changed = ('erase',t,corrected)
                if accept('erasure_correction',original[:i]+[changed]+original[i+1:],nq,
                          mutation,True):
                    phase_count += 1
                    break
        if phase_count >= 4:
            break

    # Exact equality predicates are computed BEFORE source and do not change
    # its operations. Extra flag survives only on the designated rare input.
    rare_required={case['label'] for case in rare_cases}
    rare_covered=set()
    for target in output_wires:
        for case in rare_cases:
            match = case['match']
            if not match or not set(match) <= set(ins):
                raise ValueError('Invalid rare predicate')
            literals = []
            for name,value in match.items():
                if type(value) is not int or not 0 <= value < 1 << len(ins[name]):
                    raise ValueError('Predicate outside the word width')
                literals.extend((q,(value>>i)&1) for i,q in enumerate(ins[name]))
            flips = [q for q,v in literals if not v]
            prefix = [('x',q) for q in flips]
            acc = literals[0][0]
            ladder = []
            for j,(q,_) in enumerate(literals[1:]):
                t = nq+j
                prefix.append(('and',acc,q,t))
                ladder.append((acc,q,t))
                acc = t
            flag = nq+len(ladder)
            prefix.append(('cx',acc,flag))
            prefix.extend(('erase',t,(('cz',a,b),)) for a,b,t in reversed(ladder))
            prefix.extend(('x',q) for q in reversed(flips))
            mutation=dict(label=case['label'],match=match,flag=flag,target=target)
            for key in ('semantic_roles','arithmetic_checks','arithmetic_coverage_notes'):
                if key in case:mutation[key]=case[key]
            if accept('rare_input',prefix+original+[('cx',flag,target)],flag+1,
                      mutation,True,
                      predicate=(match,flag,target)):
                rare_covered.add(case['label'])
            if len(mutants) >= max(12,(count+1)//2) and rare_covered==rare_required:
                break
        if len(mutants) >= max(12,(count+1)//2) and rare_covered==rare_required:
            break

    # General structural controls follow the targeted quota.
    for i,op in enumerate(original):
        if op[0] in ('x','cx'):
            accept('remove_gate',original[:i]+original[i+1:],nq,
                   dict(position=i,removed=op),False)
        if len(mutants)>=count:
            break
    for q in output_wires:
        if len(mutants)>=count:
            break
        accept('output_flip',original+[('x',q)],nq,dict(wire=q),False)
    for a,b in combinations(output_wires,2):
        if len(mutants)>=count:
            break
        accept('output_flip_pair',original+[('x',a),('x',b)],nq,dict(wires=[a,b]),False)
    mutants = mutants[:count]
    kept_rare={m['mutation']['label'] for m in mutants if m['kind']=='rare_input'}
    if correction_groups:
        kept_groups=[m['mutation'].get('source_segment') for m in mutants if m['kind']=='erasure_correction']
        expected_groups=[group['label'] for group in correction_groups]
        if sorted(kept_groups)!=sorted(expected_groups):
            raise ValueError('Not all different source groups of phase corrections are kept in the final quota')
    if not rare_required<=kept_rare:
        raise ValueError('Not all given rare cases are covered: '+str(sorted(rare_required-kept_rare)))
    if len(mutants)<count or sum(m['targeted'] for m in mutants)*2<count:
        raise ValueError(f'Not enough targeted canaries: {len(mutants)}, '
                         f'targeted {sum(m["targeted"] for m in mutants)}')
    return mutants
