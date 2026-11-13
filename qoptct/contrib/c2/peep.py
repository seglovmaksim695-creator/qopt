"""Пипхол: сокращение пар взаимно обратных вентилей через коммутирующие (правило Z/X-типов по кубитам).
Элементы: ('G', имя, кубиты) и ('M', анцилла, поправки) — измерение с условной поправкой (барьер по анцилле)."""
INV = {'cx': 'cx', 'cz': 'cz', 'h': 'h', 'x': 'x', 'z': 'z', 't': 'tdg', 'tdg': 't', 's': 'sdg', 'sdg': 's'}
ZT = {'z', 's', 'sdg', 't', 'tdg'}


def types(it):
    if it[0] == 'M':
        d = {it[1]: 'O'}
        for g, q in it[2]:
            for x in q:
                if x == it[1]: continue
                d[x] = 'Z' if d.get(x, 'Z') == 'Z' else 'O'
        return d
    _, g, q = it
    if g == 'cx': return {q[0]: 'Z', q[1]: 'X'}
    if g in ZT or g == 'cz': return {x: 'Z' for x in q}
    if g == 'x': return {q[0]: 'X'}
    return {x: 'O' for x in q}


def comm(a, b):
    ta, tb = types(a), types(b)
    for x in ta:
        if x in tb and (ta[x] == 'O' or tb[x] == 'O' or ta[x] != tb[x]): return False
    return True


def cancels(a, b):
    if a[0] != 'G' or b[0] != 'G': return False
    if INV.get(a[1]) != b[1]: return False
    if a[1] == 'cz': return set(a[2]) == set(b[2])
    return a[2] == b[2]


def merge_t(a, b):
    """t·t = s, tdg·tdg = sdg на одном кубите (уменьшает T)."""
    if a[0] == 'G' and b[0] == 'G' and a[2] == b[2] and a[1] == b[1] and a[1] in ('t', 'tdg'):
        return ('G', 's' if a[1] == 't' else 'sdg', a[2])
    return None


def run(items, maxlook=200):
    items = list(items); changed = True
    while changed:
        changed = False; i = 0
        while i < len(items):
            a = items[i]; j = i + 1; done = False
            while j < len(items) and j - i < maxlook:
                b = items[j]
                if cancels(a, b):
                    del items[j]; del items[i]; changed = done = True; break
                m = merge_t(a, b)
                if m:
                    del items[j]; items[i] = m; changed = done = True; break
                if not comm(a, b): break
                j += 1
            if not done: i += 1
    return items
