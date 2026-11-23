"""The full contract of the new language on top of the UNCHANGED kprove.

No wire is excluded from the comparison. Vanished LEFT registers and
freed ancillas must become zero. Junk is part of the ordinary outputs,
and the symbolic contract must specify all its bits.
"""
from hashlib import sha256
from pathlib import Path
import importlib
import json
import resource
import sys
import time


def digest_json(value):
    return sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()


def check_layout(layout,nq):
    if type(nq) is not int or nq<=0 or layout.get('nq')!=nq:
        raise ValueError('The full size of the contract does not match')
    if layout.get('garbage') or layout.get('free_wires'):
        raise ValueError('Undefined outputs are forbidden')
    if type(layout.get('partial_gate_domains_required')) is not bool:
        raise ValueError('The explicit policy of partial domains is missing')
    policy=layout.get('domain_policy')
    if policy not in ('exact_logical_AND_ERASE_v1','total_measurement_instruments_v1'):
        raise ValueError('Unknown domain policy of the source rules')
    if policy=='total_measurement_instruments_v1' and (layout['partial_gate_domains_required'] or
            layout.get('source_backend')!='qualtran_bloqgraph_measurement'):
        raise ValueError('Wrong policy of the measurement source')
    for side in ('inputs','outputs'):
        if not isinstance(layout.get(side),dict):
            raise ValueError(f'No explicit layout {side}')
        wires=[q for word in layout[side].values() for q in word]
        if any(type(q) is not int or not 0<=q<nq for q in wires) or len(set(wires))!=len(wires):
            raise ValueError(f'Overlapping or invalid wire {side}')
    if not layout['inputs'] or not layout['outputs']:
        raise ValueError('An empty signature is not supported by this client')
    free={q for word in layout['inputs'].values() for q in word}
    zeros=layout.get('initial_zero')
    if not isinstance(zeros,list) or any(type(q) is not int for q in zeros) or len(set(zeros))!=len(zeros):
        raise ValueError('Invalid set of initial zeros')
    if set(zeros)!=set(range(nq))-free:
        raise ValueError('The initial zeros are not exactly the complement of the inputs')
    return free


def prove_program(qasm,layout,symbolic,*,pre=None,validation,contract,budget=60,max_mb=900,
                  negative_control=False,domains=None):
    """symbolic(S)->{output wire: BDD}; pre(S)->BDD (a non-empty promise).

    validation receipt: status=PASS, qasm_sha256, contract_sha256, layout_sha256,
    independent_inputs>0, all_wires_checked=True, phase_checked=True.
    For negative_control the receipt has status=WITNESS, the own hash
    of the mutant, baseline_qasm_sha256 and an independent witness. Such a mutant also
    REACHES the real prover: a refusal of the external admission does not count
    as a negative control of the prover itself.
    This is a link to the external check, not a replacement of the independent check.
    The contract passed in holds a readable formula and the origin of the functions.
    For the logical policy the positive run itself checks the exact domain
    with the same pre. An external domains does not replace this check. A negative
    control is passed to kprove anyway, even if it violates the partial domain.
    """
    from qoptct.verify.kprove import kprove
    from qoptct.verify._parse import parse_qasm_text
    from qoptct.verify.ksym import SpecContext
    from qoptct.verify.bdd import ZERO
    started=time.monotonic()
    qhash=sha256(qasm.encode()).hexdigest()
    chash=digest_json(contract)
    result={}
    invoked=False
    artifact=None
    checked_domains=None
    try:
        nq,_,_=parse_qasm_text(qasm)
        free=check_layout(layout,nq)
        for key,want in dict(status='WITNESS' if negative_control else 'PASS',qasm_sha256=qhash,contract_sha256=chash,
                             layout_sha256=digest_json(layout),
                             all_wires_checked=True,phase_checked=True).items():
            if validation.get(key)!=want:
                raise ValueError(f'No matching independent check: {key}')
        if type(validation.get('independent_inputs')) is not int or validation['independent_inputs']<=0:
            raise ValueError('The independent check has no inputs')
        if negative_control:
            base=validation.get('baseline_qasm_sha256','')
            if not isinstance(base,str) or len(base)!=64 or base==qhash or not validation.get('witness'):
                raise ValueError('No baseline and no independent witness of the canary')
        if negative_control:
            # Mutants may deliberately leave the source translator's language.
            # They still reach the actual prover; the separate raw witness is
            # required above. Never count our policy gate as a prover canary.
            pass
        elif layout['domain_policy']=='exact_logical_AND_ERASE_v1':
            from .sim_parse import parse_qasm
            _,logical_ops=parse_qasm(qasm)
            required=any(op[0] in ('and','erase') for op in logical_ops)
            if layout['partial_gate_domains_required']!=required:
                raise ValueError('The partial-domain policy does not match the actual text')
            if not negative_control:
                from .domain_guard import check_domains
                import os
                checked_domains=check_domains(qasm,layout,pre=pre,contract=contract,budget=budget,max_mb=max_mb,max_nodes=int(os.environ.get('BRIDGE_DOMAIN_MAX_NODES','700000')))  # resource knob only; the default equals the 0.6 bridge
                if checked_domains['status']!='PASS':
                    raise ValueError('The exact domains of the source gates are not proved: '+checked_domains['detail'])
                if domains is not None and any(domains.get(k)!=checked_domains.get(k) for k in
                        ('status','qasm_sha256','contract_sha256','layout_sha256','pre_sha256')):
                    raise ValueError('The external domain receipt does not match the actual check')
        else:
            from .gf_measure import parse_qasm as parse_measurement_qasm
            parse_measurement_qasm(qasm)
        outwires={q for word in layout['outputs'].values() for q in word}
        clean=set(range(nq))-free
        def full_spec(S):
            out=symbolic(S)
            if not isinstance(out,dict) or set(out)!=outwires:
                raise ValueError('The specification does not cover exactly ALL source outputs and junk')
            full={q:ZERO for q in range(nq)}
            full.update(out)
            return full
        def promise(ps,xedge):
            value=pre(SpecContext(ps.B,xedge,clean))
            if value==ZERO:
                raise ValueError('An empty precondition domain is forbidden')
            if checked_domains is not None:
                from .domain_guard import bdd_input_digest
                levels={q:ps.B.level(xedge(q)) for q in free}
                if bdd_input_digest(ps.B,value,levels)!=checked_domains['pre_sha256']:
                    raise ValueError('The precondition changed between the domain guard and kprove')
            return value
        ins=layout['inputs']
        order=[word[i] for i in range(max(map(len,ins.values()))) for word in ins.values() if i<len(word)]
        from .._cache import cache_dir
        folder=cache_dir('prover_inputs')
        artifact=folder/f'{qhash}.qasm'
        if artifact.exists():
            if artifact.read_text()!=qasm:
                raise ValueError('The file with this hash holds another text')
        else:
            artifact.write_text(qasm)
        invoked=True
        result=kprove(str(artifact),spec=full_spec,n=nq,clean_wires=sorted(clean),
                      pre=promise if pre else None,order=order,engine='bdd',
                      sample=False,confirm=False,budget_sec=budget,max_mb=max_mb)
        result['prover_invoked']=True
    except Exception as exc:
        result=dict(status='ERROR',detail=f'{type(exc).__name__}: {exc}',prover_invoked=invoked)
    module=importlib.import_module('qoptct.verify.kprove')
    path=Path(module.__file__).resolve()
    # kprove may return internal callables for optional witness confirmation.
    result={k:v for k,v in result.items() if not callable(v)}
    result.update(qasm_sha256=qhash,contract_sha256=chash,
                  prover_file='qoptct/verify/kprove.py',prover_sha256=sha256(path.read_bytes()).hexdigest(),
                  wall_seconds=time.monotonic()-started,
                  peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1e6 if sys.platform=='darwin' else 1000),
                  method='bdd',complete_output_contract=True)
    if checked_domains is not None:
        result['partial_domains']=checked_domains
    if artifact is not None:
        result['artifact']='prover_inputs/'+artifact.name
    return result
