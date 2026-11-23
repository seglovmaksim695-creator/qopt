"""Exact domains of the partial AND/ERASE gates, independent of the final equivalence.

A client of the unchanged ROBDD. Never returns an arithmetic PROVED. PASS means
that at every actual step the AND target is 0 and the erased bit equals
the given correction EXACTLY, constant included (not only up to a phase).
"""
from hashlib import sha256
import time

from .sim_parse import parse_qasm
from .strict_proof import check_layout,digest_json


def bdd_input_digest(B,edge,wire_levels):
    """Canonical function with the names of the physical inputs, not node numbers."""
    from qoptct.verify.bdd import ONE,ZERO
    by_level={level:wire for wire,level in wire_levels.items()}
    memo={ONE:'1',ZERO:'0'}
    def visit(e):
        if e not in memo:
            level=B.level(e)
            if level not in by_level:
                raise ValueError('The precondition uses a variable that is not an input')
            hi,lo=B.cof(e,level)
            memo[e]=sha256(f'q{by_level[level]}:{visit(hi)}:{visit(lo)}'.encode()).hexdigest()
        return memo[e]
    return sha256(visit(edge).encode()).hexdigest()


def check_domains(qasm,layout,*,pre=None,contract,budget=60,max_nodes=700000,max_mb=800):
    from qoptct.verify.bdd import BDD,ONE,ZERO
    from qoptct.verify.ksym import SpecContext
    from qoptct.verify.kprove import MemGuard
    started=time.monotonic()
    result=dict(status='UNKNOWN',scope='exact partial-gate domains, not arithmetic',
                qasm_sha256=sha256(qasm.encode()).hexdigest(),contract_sha256=digest_json(contract),
                layout_sha256=digest_json(layout))
    checks=0
    try:
        nq,ops=parse_qasm(qasm)
        free=check_layout(layout,nq)
        B=BDD()
        mem=MemGuard(max_mb)
        initial=[ZERO]*nq
        variables={}
        ins=layout['inputs']
        for i in range(max(map(len,ins.values()))):
            for word in ins.values():
                if i<len(word):
                    level,edge=B.new_var()
                    variables[word[i]]=level
                    initial[word[i]]=edge
        def guard():
            if time.monotonic()-started>budget or B.nodes()>max_nodes:
                raise TimeoutError('Limit of the exact domain check')
            mem.check()
        B.hook=guard
        P=pre(SpecContext(B,lambda q:initial[q],set(range(nq))-free)) if pre else ONE
        if P==ZERO:
            raise ValueError('Empty promised domain')
        result['pre_sha256']=bdd_input_digest(B,P,variables)
        state=initial.copy()
        for i,op in enumerate(ops):
            if i%128==0:
                guard()
            kind=op[0]
            bad=ZERO
            if kind=='x':
                state[op[1]] ^= 1
            elif kind=='cx':
                state[op[2]]=B.XOR(state[op[2]],state[op[1]])
            elif kind=='and':
                _,a,b,t=op
                bad=B.AND(P,state[t])
                state[t]=B.AND(state[a],state[b])
                checks+=1
            elif kind=='erase':
                _,t,corrections=op
                phase=state[t]
                for correction in corrections:
                    if correction[0]=='z':
                        phase=B.XOR(phase,state[correction[1]])
                    else:
                        phase=B.XOR(phase,B.AND(state[correction[1]],state[correction[2]]))
                bad=B.AND(P,phase)
                state[t]=ZERO
                checks+=1
            else:
                raise ValueError(f'Unsupported exact domain atom: {kind}')
            if bad!=ZERO:
                assignment=B.pick_sat(bad)
                values={name:sum(int(assignment.get(variables[q],0))<<j for j,q in enumerate(word))
                        for name,word in ins.items()}
                result.update(status='REFUTED',operation_index=i,operation=op,
                              witness=values,detail='The exact domain of a partial gate is violated')
                break
        else:
            result.update(status='PASS',detail='All AND/ERASE have their exact domain at every step')
        result.update(bdd_nodes=B.nodes(),domain_checks=checks)
    except Exception as exc:
        result.update(status='UNKNOWN' if isinstance(exc,TimeoutError) or type(exc).__name__=='Budget' else 'ERROR',
                      detail=f'{type(exc).__name__}: {exc}')
    result['wall_seconds']=time.monotonic()-started
    return result
