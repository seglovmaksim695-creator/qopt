# SPDX-License-Identifier: Apache-2.0
"""E-transform: replace the final CX block N (followed by a diagonal/measurement tail) by a cheaper N' that is exact on data rows and
changes garbage (X-measured) rows by E = L L'^-1 (unit upper triangular in measurement order + data columns); rewrite the tail accordingly."""
import random
import time

from .ir import pcond, fcond, counts

def bm(cxs,n):
    S=[1<<i for i in range(n)]
    for c,t in cxs: S[t]^=S[c]
    return S
def inv_rows(rows,n):
    A=list(rows); I=[1<<k for k in range(n)]
    for c in range(n):
        p=next((r for r in range(c,n) if A[r]>>c&1),None)
        if p is None: return None
        A[c],A[p]=A[p],A[c]; I[c],I[p]=I[p],I[c]
        for r in range(n):
            if r!=c and A[r]>>c&1: A[r]^=A[c]; I[r]^=I[c]
    return I
def mul(Aa,Bb):
    out=[]
    for r in Aa:
        v=0;k=0
        while r:
            if r&1: v^=Bb[k]
            r>>=1;k+=1
        out.append(v)
    return out
def locate(ops):
    fm=next(i for i,o in enumerate(ops) if o[0]=='M')
    i=fm-1
    while i>=0 and not (ops[i][0]=='G' and ops[i][1]=='cx'): i-=1
    bend=i; j=i
    while j>=0 and ops[j][0]=='G' and ops[j][1]=='cx': j-=1
    return j,bend   # block = ops[j+1:bend+1]
class Bail(Exception): pass
def analyse(ops,n):
    fm=next((i for i,o in enumerate(ops) if o[0]=='M'),None)
    if fm is None or fm==0: raise Bail('no measurement')
    hp=fm-1
    if not (ops[hp][0]=='G' and ops[hp][1]=='h'): raise Bail('no X measurement')
    # tail starts at first op after region; region = maximal run before the (x*) h of first M
    i=hp-1
    while i>=0 and ops[i][0]=='G' and ops[i][1] in('x',): i-=1
    # i is last op of region candidates
    e=i
    while i>=0 and ops[i][0]=='G' and ops[i][1] in('cx','s','sdg','z','cz'): i-=1
    j=i; bend=e
    blk=[];regterms=[];S0=[1<<k for k in range(n)]
    for o in ops[j+1:bend+1]:
        nm,qs=o[1],o[2]
        if nm=='cx': S0[qs[1]]^=S0[qs[0]]; blk.append(qs)
        else: regterms.append((nm,qs,tuple(S0[q] for q in qs)))
    if not blk: raise Bail('no block')
    tail=ops[bend+1:]
    meas=[];bits={}
    for t,o in enumerate(tail):
        if o[0]=='R': raise Bail('reset in tail')
        if o[0]=='G' and o[1] not in('x','s','sdg','z','cz','h'): raise Bail('gate '+o[1])
        if o[0]=='G' and o[1]=='h':
            if not(t+1<len(tail) and tail[t+1][0]=='M' and tail[t+1][1]==o[2][0]): raise Bail('h without M')
        if o[0]=='M': meas.append(o[1]); bits[o[1]]=o[2]
        if o[0]=='IF':
            for g,qs in o[2]:
                if g not in('cz','s','sdg','z','x'): raise Bail('IF gate '+g)
    if len(set(meas))!=len(meas): raise Bail('wire measured twice')
    return j,bend,blk,tail,meas,bits,regterms
def allowed_masks(data,meas,n):
    order={w:k for k,w in enumerate(meas)}
    allowed={}
    for w in meas:
        m=0
        for d in data: m|=1<<d
        for l in meas:
            if order[l]>=order[w]: m|=1<<l
        allowed[w]=m
    return allowed
def valid(L,S,n,data,meas,allowed):
    if any(L[d]!=S[d] for d in data): return False
    Li=inv_rows(L,n)
    if Li is None: return False
    E=mul(S,Li)
    return all(E[w]&~allowed[w]==0 for w in meas)
def delete_greedy(blk,n,data,meas,rng=None,passes=50):
    S=bm(blk,n); allowed=allowed_masks(data,meas,n)
    cur=list(blk)
    changed=True
    while changed:
        changed=False
        idx=list(range(len(cur)-1,-1,-1))
        if rng: rng.shuffle(idx)
        for k in sorted(idx,reverse=True) if not rng else idx:
            if k>=len(cur): continue
            t=cur[:k]+cur[k+1:]
            if valid(bm(t,n),S,n,data,meas,allowed): cur=t; changed=True
    return cur
def rewrite(ops,n,j,bend,blk,Np,tail,meas,bits,rt=()):
    data=[w for w in range(n) if w not in meas]; G=len(meas)
    L=bm(blk,n); Lp=bm(Np,n); Li=inv_rows(Lp,n); E=mul(L,Li)
    gidx={w:k for k,w in enumerate(meas)}
    # A and T
    A=[[(E[meas[a]]>>meas[b])&1 for b in range(G)] for a in range(G)]
    # invert A (G x G) -> A^{-1}; T = (A^-1)^T
    Arows=[sum(A[a][b]<<b for b in range(G)) for a in range(G)]
    Ai=inv_rows(Arows,G)
    T=[[ (Ai[l]>>jj)&1 for l in range(G)] for jj in range(G)]   # T[j][l] = (A^-1)^T[j][l] = Ai[l][j]
    for jj in range(G):
        for l in range(G):
            if T[jj][l] and l>jj: raise Bail('T not lower triangular')
    bitlist=[bits[w] for w in meas]
    def sub(bset):
        res=set()
        for b in bset:
            if b=='1': res^={'1'}; continue
            jj=bitlist.index(b)
            res^={bitlist[l] for l in range(G) if T[jj][l]}
        return frozenset(res)
    def cstr(s):
        s=sorted(x for x in s if x!='1')
        return '('+' ^ '.join(f'c[{b}]' for b in s)+')'
    # constants kappa
    kappa=[0]*n; first=0
    while first<len(tail) and tail[first][0]=='G' and tail[first][1]=='x':
        kappa[tail[first][2][0]]^=1; first+=1
    expr={}
    for w in data: expr[w]=(1<<w,0)
    for w in meas:
        c=0
        for k in range(n):
            if E[w]>>k&1 and k in data and kappa[k]: c^=1
        c^=kappa[w]
        expr[w]=(E[w],c)
    def supp(m): return [k for k in range(n) if m>>k&1]
    Li_x=Li
    out=list(ops[:j+1])+[('G','cx',p) for p in Np]
    for t in range(first): out.append(tail[t])
    pend={}   # (k,l)-> frozenset condset
    measured=set()
    def xorset(a,b): return frozenset(set(a)^set(b))
    def emit_cz_cond(k,l,cs):
        if '1' in cs: out.append(('G','cz',(k,l)))
        rest=[x for x in cs if x!='1']
        if rest: out.append(('IF',cstr(rest),[('cz',(k,l))]))
    def flush(wire=None):
        for key in list(pend):
            if wire is None or wire in key:
                cs=pend.pop(key)
                if cs: emit_cz_cond(key[0],key[1],cs)
    def add_pair(k,l,cs):
        key=(min(k,l),max(k,l))
        pend[key]=xorset(pend.get(key,frozenset()),cs)
    def emit_lin(gate,k,groups):
        # groups: list of frozensets (new bits) AND; '1' means unconditional (skip)
        gs=[g for g in groups if g!=frozenset({'1'})]
        if any(len(g)==0 for g in gs): return
        if not gs: out.append(('G',gate,(k,)))
        else: out.append(('IF',' && '.join(cstr(g) for g in gs),[(gate,(k,))]))
    xend=[]
    for nm,qs,snaps in rt:
        ex=[]
        for m in snaps:
            mm=0
            for k in range(n):
                if m>>k&1: mm^=Li_x[k]
            ex.append((mm,0))
        gsub1=[frozenset({'1'})]
        if nm=='z':
            for k in supp(ex[0][0]): emit_lin('z',k,gsub1)
        elif nm in('s','sdg'):
            sg=1 if nm=='s' else -1; S_=supp(ex[0][0])
            for k in S_: emit_lin('s' if sg>0 else 'sdg',k,gsub1)
            for x_ in range(len(S_)):
                for y_ in range(x_+1,len(S_)): add_pair(S_[x_],S_[y_],gsub1[0])
        elif nm=='cz':
            Sa,Sb=supp(ex[0][0]),supp(ex[1][0]); cnt={}
            for k in Sa:
                for l in Sb:
                    key=(k,l) if k<l else (l,k); cnt[key]=cnt.get(key,0)^1
            for (k,l),v in cnt.items():
                if not v: continue
                if k==l: emit_lin('z',k,gsub1)
                else: add_pair(k,l,gsub1[0])
    tpos=first
    while tpos<len(tail):
        o=tail[tpos]; tpos+=1
        if o[0]=='G' and o[1]=='h':
            w=o[2][0]; m=tail[tpos]; tpos+=1
            for k in supp(expr[w][0]):
                if k in measured and k!=w: pass
            flush(w); out.append(o); out.append(m); measured.add(w); continue
        if o[0]=='IF' and all(g=='x' for g,_ in o[2]):
            gg=[sub(frozenset(g)) for g in pcond(o[1])]
            if any(not g for g in gg): continue
            xend.append(('IF',' && '.join(cstr(g) for g in gg),o[2])); continue
        if o[0]=='G': groups=[frozenset({'1'})]; body=[(o[1],o[2])]
        else:
            groups=[frozenset(g) for g in pcond(o[1])]; body=o[2]
        gsub=[sub(g) for g in groups]
        for g,qs in body:
            if g=='x': raise Bail('x in mixed IF')
            for q in qs:
                if q in meas and q in measured: raise Bail('term on measured wire')
            if g=='z':
                for k in supp(expr[qs[0]][0]): emit_lin('z',k,gsub)
            elif g in('s','sdg'):
                if len(gsub)!=1: raise Bail('S with AND cond')
                mk,kp=expr[qs[0]]; sg=(1 if g=='s' else -1)*(-1 if kp else 1)
                S_=supp(mk)
                for k in S_: emit_lin('s' if sg>0 else 'sdg',k,gsub)
                for x in range(len(S_)):
                    for y in range(x+1,len(S_)): add_pair(S_[x],S_[y],gsub[0])
            elif g=='cz':
                if len(gsub)!=1: raise Bail('cz with AND cond')
                (ma,ka),(mb,kb)=expr[qs[0]],expr[qs[1]]
                Sa,Sb=supp(ma),supp(mb)
                cnt={}
                for k in Sa:
                    for l in Sb:
                        key=(k,l) if k<l else (l,k)
                        cnt[key]=cnt.get(key,0)^1
                for (k,l),v in cnt.items():
                    if not v: continue
                    if k==l: emit_lin('z',k,gsub)
                    else: add_pair(k,l,gsub[0])
                if ka:
                    for l in Sb: emit_lin('z',l,gsub)
                if kb:
                    for k in Sa: emit_lin('z',k,gsub)
    flush()
    # B part: data wire k gets Z with cond XOR of original bits of garbage j with E[w_j] bit k
    for k in data:
        cs=set()
        for jj,w in enumerate(meas):
            if E[w]>>k&1: cs^={bits[w]}
        cs=sub(cs)
        if cs: out.append(('IF',cstr(cs),[('z',(k,))]))
    out+=xend
    return out
def cost_of(ops,n,j,bend,blk,Np,tail,meas,bits,rt=()):
    try:
        out=rewrite(ops,n,j,bend,blk,Np,tail,meas,bits,rt)
    except Bail:
        return 10**9,None
    return counts(out)[1],out
def optimize(ops,n,seed=0,rounds=6,verbose=False,tlimit=1e9):
    """Greedy deletion search for N' by total CX (block + rewritten corrections). Returns the best ops list."""
    t0=time.time()
    j,bend,blk,tail,meas,bits,rt=analyse(ops,n)
    data=[w for w in range(n) if w not in meas]
    S=bm(blk,n); allowed=allowed_masks(data,meas,n)
    rng=random.Random(seed)
    cur=list(blk); c0,out=cost_of(ops,n,j,bend,blk,cur,tail,meas,bits,rt)

    best=(c0,list(cur),out)
    for r in range(rounds):
        idx=list(range(len(cur))); rng.shuffle(idx); changed=False
        for k in sorted(idx,reverse=True) if r==0 else idx:
            if time.time()-t0>tlimit: break
            if k>=len(cur): continue
            t=cur[:k]+cur[k+1:]
            if not valid(bm(t,n),S,n,data,meas,allowed): continue
            c,o=cost_of(ops,n,j,bend,blk,t,tail,meas,bits,rt)
            if c<=c0:
                if c<c0: changed=True
                cur=t; c0=c; out=o
                if c<best[0]: best=(c,list(cur),o)
        if verbose: print('round',r,'N\'',len(cur),'total',c0,flush=True)
        if not changed and r>0: break
    if best[2] is None: raise Bail('no valid rewrite')
    return best[2]
