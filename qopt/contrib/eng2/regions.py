"""Разбиение схемы на области без H: {cx,x,s,sdg,t,tdg}. Область ставится на месте закрытия
(все её кубиты к тому моменту не затронуты чужими вентилями), вентили вне области коммутируют с ней (непересекающиеся кубиты)."""
NOH = {"cx","x","s","sdg","t","tdg"}
def split_regions(gates, n):
    cur=[None]*n; regs={}; out=[]; nid=0
    def close(r):
        out.append(('R',r, sorted(regs[r]['q']), [g for _,g in sorted(regs[r]['g'])]))
        for q in regs[r]['q']: cur[q]=None
    for i,g in enumerate(gates):
        nm,qs=g
        if nm not in NOH:
            for q in qs:
                if cur[q] is not None: close(cur[q])
            out.append(('G',g)); continue
        rs={cur[q] for q in qs if cur[q] is not None}
        if not rs:
            r=nid; nid+=1; regs[r]={'q':set(),'g':[]}
        else:
            r=min(rs)
            for o in rs:
                if o!=r:
                    regs[r]['q']|=regs[o]['q']; regs[r]['g']+=regs[o]['g']
        regs[r]['g'].append((i,g)); regs[r]['q']|=set(qs)
        for q in regs[r]['q']: cur[q]=r
    for r in sorted({c for c in cur if c is not None}): close(r)
    return out


def split_regions2(gates, n):
    """Области с поэтапным закрытием кубитов: H на q закрывает только q; область «запечатывается» (ставится в позицию первого H),
    после этого принимает только вентили на своих ещё открытых кубитах."""
    cur = [None] * n; regs = {}; out = []; nid = 0
    def seal(r):
        if not regs[r]['sealed']:
            regs[r]['sealed'] = True; out.append(('P', r))
    for i, g in enumerate(gates):
        nm, qs = g
        if nm not in NOH:  # h
            q = qs[0]; r = cur[q]
            if r is not None:
                seal(r); cur[q] = None
            out.append(('G', g)); continue
        touched = {cur[q] for q in qs if cur[q] is not None}
        sealed = [r for r in touched if regs[r]['sealed']]
        if sealed:
            if len(touched) == 1 and all(cur[q] == sealed[0] for q in qs):
                regs[sealed[0]]['g'].append((i, g)); continue
            for q in qs:
                if cur[q] is not None and regs[cur[q]]['sealed']: cur[q] = None
            touched = {cur[q] for q in qs if cur[q] is not None}
        if not touched:
            r = nid; nid += 1; regs[r] = {'q': set(), 'g': [], 'sealed': False}
        else:
            r = min(touched)
            for o in touched:
                if o != r:
                    regs[r]['q'] |= regs[o]['q']; regs[r]['g'] += regs[o]['g']; regs[o]['g'] = []; regs[o]['dead'] = True
        regs[r]['g'].append((i, g)); regs[r]['q'] |= set(qs)
        for q in regs[r]['q']: cur[q] = r
    for r in sorted({c for c in cur if c is not None}): seal(r)
    res = []
    for it in out:
        if it[0] == 'G': res.append(it)
        else:
            R = regs[it[1]]
            if R.get('dead'): continue
            gl = [g for _, g in sorted(R['g'])]
            res.append(('R', it[1], sorted({q for _, qq in gl for q in qq}), gl))
    return res
