"""Canonical text. No circuit recognition and no reordering."""
AND10 = (('h',2),('t',2),('cx',1,2),('tdg',2),('cx',0,2),
         ('t',2),('cx',1,2),('tdg',2),('h',2),('sdg',2))

def emit_qasm(nq, ops):
    lines = ['OPENQASM 3;', 'include "stdgates.inc";', f'qubit[{nq}] q;',
             f'bit[{max(1,sum(op[0]=="erase" for op in ops))}] c;']
    k = 0
    def gate(g, wires):
        if len(set(wires)) != len(wires) or any(type(w) is not int or not 0 <= w < nq for w in wires):
            raise ValueError(f'Invalid wires {g} {wires}')
        return g + ' ' + ', '.join(f'q[{w}]' for w in wires) + ';'
    for op in ops:
        if op[0] in ('x','cx'):
            lines.append(gate(op[0],op[1:]))
        elif op[0]=='and':
            a,b,t=op[1:]
            gate('and',(a,b,t))
            for g,*w in AND10:
                lines.append(gate(g,[(a,b,t)[i] for i in w]))
        elif op[0]=='erase':
            _,t,corr=op
            lines += [gate('h',(t,)),f'c[{k}] = measure q[{t}];']
            body=[]
            for g,*w in corr:
                if g not in ('z','cz') or t in w:
                    raise ValueError('Unsupported erase correction')
                body.append(gate(g,w))
            body.append(gate('x',(t,)))
            lines.append(f'if (c[{k}]) '+'{'+' '.join(body)+'}')
            k+=1
        else:
            raise ValueError(f'Unknown operation {op}')
    return '\n'.join(lines)+'\n'
