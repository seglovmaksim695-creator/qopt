"""Copy of bridge/rules.py with one addition: an opaque leaf -> the marker ('opaque', leaf index).
The other rules are byte for byte the same."""
"""Every leaf gets exactly one explicit rule; the history is kept."""
from qoptct.foreign.qualtran.bridge.emit import emit_qasm

from qoptct.foreign.qualtran.bridge.rules import TranslationError

def translate(lowered, inputs, outputs):
    nq=lowered['nq']
    if 'inputs' in lowered and inputs!=lowered['inputs']:
        raise TranslationError('The original input signature of Qualtran must not be narrowed')
    if 'outputs' in lowered and outputs!=lowered['outputs']:
        raise TranslationError('The original output signature of Qualtran must not be changed')
    for kind,regs in (('inputs',inputs),('outputs',outputs)):
        flat=[q for w in regs.values() for q in w]
        if len(flat)!=len(set(flat)) or any(type(q) is not int or not 0<=q<nq for q in flat):
            raise TranslationError(f'Overlapping or invalid layout {kind}')
    free={q for w in inputs.values() for q in w}
    if any(not 0<=q<nq for q in free):
        raise TranslationError('Input outside the original circuit')
    constants=[None if q in free else 0 for q in range(nq)]
    ops=[]; trace=[]; scratch=None
    def x(t):
        ops.append(('x',t))
        if constants[t] is not None: constants[t]^=1
    def cx(a,t):
        ops.append(('cx',a,t))
        if constants[a]==0: return
        if constants[a]==1 and constants[t] is not None: constants[t]^=1
        else: constants[t]=None
    def AND(a,b,t):
        if constants[t]!=0:
            raise TranslationError(f'The purity of the target is not set for And q[{t}]')
        ops.append(('and',a,b,t))
        constants[t]=0 if 0 in (constants[a],constants[b]) else (1 if constants[a]==constants[b]==1 else None)
    def erase(a,b,t):
        ops.append(('erase',t,(('cz',a,b),)))
        constants[t]=0
    for i,leaf in enumerate(lowered['leaves']):
        g=leaf['name']; qs=leaf['qubits']; start=len(ops)
        obligations=[]
        try:
            if any(type(q) is not int or not 0<=q<nq for q in qs) or len(set(qs))!=len(qs):
                raise TranslationError('Repeated or invalid wire')
            if g=='x' and len(qs)==1: x(*qs)
            elif g=='cx' and len(qs)==2: cx(*qs)
            elif g=='swap' and len(qs)==2:
                a,b=qs; old=constants[a],constants[b]
                cx(a,b);cx(b,a);cx(a,b)
                constants[a],constants[b]=old[1],old[0]
            elif g=='ccx' and len(qs)==3:
                if scratch is None:
                    scratch=len(constants);constants.append(0)
                a,b,t=qs
                AND(a,b,scratch);cx(scratch,t);erase(a,b,scratch)
            elif g in ('controlled_and','controlled_erase') and len(qs)==4:
                from qoptct.foreign.qualtran.bridge.controlled_and import expand
                if scratch is None:
                    scratch=len(constants)
                expanded=expand(leaf,ancilla=scratch,nq=len(constants),
                                known_zero=[q for q,v in enumerate(constants) if v==0])
                constants.extend([0]*(expanded['nq']-len(constants)))
                for primitive in expanded['ops']:
                    if primitive[0]=='x':x(primitive[1])
                    elif primitive[0]=='cx':cx(*primitive[1:])
                    elif primitive[0]=='and':AND(*primitive[1:])
                    elif primitive[0]=='erase':
                        target,corrections=primitive[1:]
                        if len(corrections)!=1 or corrections[0][0]!='cz':
                            raise TranslationError('Unexpected correction of the A1 extension')
                        erase(*corrections[0][1:],target)
                    else:raise TranslationError('Unknown atom of the A1 extension')
                obligations=expanded['domain_obligations']
            elif g in ('and','erase') and len(qs)==3:
                a,b,t=qs; cvs=leaf.get('controls',[1,1])
                if len(cvs)!=2 or any(type(v) is not int or v not in (0,1) for v in cvs):
                    raise TranslationError('Unknown control convention')
                flips=[q for q,v in zip((a,b),cvs) if v==0]
                for q in flips:x(q)
                if g=='and':AND(a,b,t)
                else:erase(a,b,t)
                for q in reversed(flips):x(q)
            elif g=='opaque':
                # ew: black box (multiplier/Kaliski) is carried over as a marker; the right outputs are unknown,
                # the LEFT wires it releases are zero after the call (part of the contract of the box).
                ops.append(('opaque',i))
                for q in qs: constants[q]=None
                for r,ws in leaf['regs_in'].items():
                    if r not in leaf['regs_out']:
                        for q in ws: constants[q]=0
            else: raise TranslationError(f'No rule for {g}')
        except Exception as e:
            raise TranslationError(f'Leaf {i}, {leaf.get("path")}, {g}{qs}: {e}') from e
        trace.append(dict(leaf=leaf,source_index=i,logical_start=start,logical_end=len(ops),
                          domain_obligations=obligations))
    N=len(constants)
    layout=dict(nq=N,source_nq=nq,inputs=inputs,outputs=outputs,
                initial_zero=[q for q in range(N) if q not in free],
                garbage=[],bridge_scratch=[] if scratch is None else [scratch],
                source_wire_names=lowered.get('wire_names',[]), bit_order='LSB-first in every layout list')
    layout['partial_gate_domains_required']=any(op[0] in ('and','erase') for op in ops)
    layout['domain_policy']='exact_logical_AND_ERASE_v1'
    for key in ('registers','source_backend','thru_permutations','retired_wires'):
        if key in lowered:
            layout[key]=lowered[key]
    return dict(nq=N,ops=ops,layout=layout,trace=trace,qasm=emit_qasm(N,ops),
                requires_exact_partial_domains=layout['partial_gate_domains_required'])


def translate_ew(lowered):
    """Like translate, but no QASM is built (it has no opaque markers)."""
    global emit_qasm
    saved = emit_qasm
    emit_qasm = lambda nq, ops: None
    try:
        return translate(lowered, lowered['inputs'], lowered['outputs'])
    finally:
        emit_qasm = saved
