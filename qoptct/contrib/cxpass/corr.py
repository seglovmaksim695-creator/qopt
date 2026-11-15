# SPDX-License-Identifier: Apache-2.0
"""Relocation of conditional CZ corrections: cz(Wi,Wj) = S(Wi) S(Wj) Sdg(Wi^Wj) (exact, on wire values).
Components are placed (as conditional free S/Sdg gates) at any later time when some data wire holds the needed form."""
def track(ops,n):
    """Intervals (start, end, wire, form, const) during which a wire holds an affine form; anc = wires that are ever
    re-initialised (H, measurement, reset, conditional non-diagonal gates) and are not used as relocation targets."""
    form=[1<<i for i in range(n)]; const=[0]*n; nl=n
    start=[-1]*n
    ivs=[]; snaps=[]
    anc=set()
    for o in ops:
        if o[0]=='G' and o[1]=='h': anc.add(o[2][0])
        if o[0] in ('M','R'): anc.add(o[1])
        if o[0]=='G' and o[1] not in ('cx','x','h','z','s','sdg','t','tdg','cz'): anc.update(o[2])
        if o[0]=='IF':
            for nm,qs in o[2]:
                if nm not in ('z','s','sdg','t','tdg','cz'): anc.update(qs)
    def close(w,i): ivs.append((start[w],i,w,form[w],const[w]))
    for i,o in enumerate(ops):
        gs=None
        if o[0]=='G':
            nm,qs=o[1],o[2]
            if nm=='cx':
                c,t=qs; close(t,i-1); form[t]^=form[c]; const[t]^=const[c]; start[t]=i
            elif nm=='x':
                t=qs[0]; close(t,i-1); const[t]^=1; start[t]=i
            elif nm=='h':
                t=qs[0]; close(t,i-1); form[t]=1<<nl; nl+=1; const[t]=0; start[t]=i
        elif o[0] in ('M','R'):
            t=o[1]; close(t,i-1); form[t]=1<<nl; nl+=1; const[t]=0; start[t]=i
        else:
            for nm,qs in o[2]:
                if nm=='x': t=qs[0]; close(t,i-1); form[t]=0; const[t]=0; start[t]=i
    for w in range(n): close(w,len(ops)-1)
    return ivs,anc
def forms_at(ops,n,idx):
    form=[1<<i for i in range(n)]; const=[0]*n; nl=n
    for i,o in enumerate(ops[:idx]):
        if o[0]=='G':
            nm,qs=o[1],o[2]
            if nm=='cx': form[qs[1]]^=form[qs[0]]; const[qs[1]]^=const[qs[0]]
            elif nm=='x': const[qs[0]]^=1
            elif nm=='h': form[qs[0]]=1<<nl; nl+=1; const[qs[0]]=0
        elif o[0] in ('M','R'): form[o[1]]=1<<nl; nl+=1; const[o[1]]=0
        else:
            for nm,qs in o[2]:
                if nm=='x': form[qs[0]]=0; const[qs[0]]=0
    return form,const
def relocate(ops,n):
    ivs,anc=track(ops,n)
    new_after={}  # idx -> list of ops to insert after op idx
    drop=set(); removed=0
    # iterate IF cz in order
    for i,o in enumerate(ops):
        if o[0]!='IF': continue
        czs=[(nm,qs) for nm,qs in o[2] if nm=='cz']
        if not czs: continue
        newbody=[(nm,qs) for nm,qs in o[2] if nm!='cz']
        form,const=forms_at(ops,n,i)
        okall=True; plans=[]
        for _,(a,b) in czs:
            fa,ka,fb,kb=form[a],const[a],form[b],const[b]
            if a in anc or b in anc: okall=False; break
            # S on Wa, S on Wb at i (now); Sdg on Wa^Wb: find interval after i
            target=fa^fb; tk=ka^kb
            found=None
            for (s,e,w,f,k) in ivs:
                if w in anc or f!=target: continue
                if e<=i: continue
                p=max(s,i)   # insert after op p (p>=i)
                if found is None or p<found[0]: found=(p,w,k)
            if found is None: okall=False; break
            p,w,k=found
            flip = (k^tk)  # wire value = target^k ; need value of Wa^Wb=target^tk
            plans.append((a,b,p,w,flip))
        if not okall: continue
        # apply plan
        body=newbody+[('s',(a,)) for a,b,p,w,f in plans]+[('s',(b,)) for a,b,p,w,f in plans]
        ops[i]=('IF',o[1],body)
        for a,b,p,w,f in plans:
            new_after.setdefault(p,[]).append(('IF',o[1],[('sdg' if not f else 's',(w,))]))
        removed+=len(plans)
    out=[]
    for i,o in enumerate(ops):
        out.append(o)
        out+=new_after.get(i,[])
    return out,removed
