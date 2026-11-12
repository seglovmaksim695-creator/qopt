"""И4б: сокращение по ВСЕЙ схеме (как проход cancel у Nam и др.): команда ищет вперёд, перепрыгивая коммутирующие,
партнёра — обратную (сократить пару) или диагональную фазу на том же кубите (слить). Без ограничения на число кубитов."""
from block_synth import commutes, DIAG
PH = {"t": 1, "s": 2, "sdg": 6, "tdg": 7}
SEQ = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}
SELF = {"cx", "x", "h"}


def one_sweep(c):
    c = list(c); i = 0; changed = False
    while i < len(c):
        g = c[i]; j = i + 1; done = False
        while j < len(c):
            h = c[j]
            if not (set(g[1]) & set(h[1])):
                j += 1; continue
            if g[0] in SELF and h == g:  # g ... g (через коммутирующие) = тождество
                del c[j]; del c[i]; changed = done = True; break
            if g[0] in DIAG and h[0] in DIAG and g[1] == h[1]:
                ph = (PH[g[0]] + PH[h[0]]) % 8
                new = [(n, g[1]) for n in SEQ[ph]]
                if len(new) < 2:
                    del c[j]; c[i:i + 1] = new; changed = done = True; break
            if commutes(g, h):
                j += 1; continue
            break
        if not done:
            i += 1
    return c, changed


def global_cancel(c):
    changed = True
    while changed:
        c, changed = one_sweep(c)
    return c
