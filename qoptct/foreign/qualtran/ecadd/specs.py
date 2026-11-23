"""Спецификации шагов ECAdd (формулы написаны по исходнику Qualtran, НЕ по схеме).

Слова — списки рёбер ROBDD, LSB первым; флаги — списки из одного ребра. O[c] — выходы c-го ПРЯМОГО ящика:
MM -> target, qrom_indices, reduced; KINV -> x (результат на месте), junk.
patched=True — исправленная версия (патч E05), False — исходная Qualtran 0.7.0.
"""
from qoptct.verify.bdd import ONE, ZERO
from .bv import (const, xor, mask, ite, lt_const, eq, is_zero, all_of, modadd, modsub, modneg, moddbl, add, lt, _s)


class _Base:
    range_words = ('a', 'b', 'x', 'y', 'lam', 'lam_r')   # диапазоны, известные до разбиения (подмножество pre)

    def __init__(self, n, p, patched=True, observe=False):
        self.n, self.p, self.patched, self.observe_on = n, p, patched, observe

    def rng(self, B, *words):
        """Список частей предусловия (по одной на слово) — для проверок по связным частям."""
        return [lt_const(B, w, self.p) for w in words]

    def calls(self, B, W, O):
        return []


class Step1(_Base):
    name = 'step1'

    def pre(self, B, W, O):
        return self.rng(B, W['a'], W['b'], W['x'], W['y'])

    def post(self, B, W, O):
        a, b, x, y = W['a'], W['b'], W['x'], W['y']
        e1 = eq(B, a, x)
        e2 = eq(B, b, modneg(B, y, self.p))
        e3 = B.AND(is_zero(B, a), is_zero(B, b))
        e4 = B.AND(is_zero(B, x), is_zero(B, y))
        if self.patched:
            ctrl = all_of(B, e3 ^ 1, e4 ^ 1, B.AND(e1, e2) ^ 1)
            f2 = B.AND(e2, ctrl ^ 1)
        else:
            ctrl = all_of(B, e2 ^ 1, e3 ^ 1, e4 ^ 1)
            f2 = e2
        return dict(f1=[e1], f2=[f2], f3=[e3], f4=[e4], ctrl=[ctrl], a=a, b=b, x=x, y=y)


class Step2(_Base):
    name = 'step2'

    def pre(self, B, W, O):
        return self.rng(B, W['a'], W['b'], W['x'], W['y'], W['lam_r'])

    def _xy(self, B, W):
        p, ctrl = self.p, W['ctrl'][0]
        x2 = modsub(B, W['x'], W['a'], p)
        y2 = ite(B, ctrl, modsub(B, W['y'], W['b'], p), W['y'])
        return x2, y2

    def calls(self, B, W, O):
        x2, y2 = self._xy(B, W)
        return [('KINV', {'x': x2}), ('MM', {'x': O[0]['x'], 'y': y2})]

    def splits(self, B, W, O):
        """Разбиение случаев (для памяти при n=256): исходы сравнений внутри ModSub/CModSub."""
        out = [('x<a', lt(B, W['x'], W['a']))]
        if W['ctrl'][0] != ZERO:
            out.append(('y<b', lt(B, W['y'], W['b'])))
        return out

    def post(self, B, W, O):
        f1, ctrl, lam_r = W['f1'][0], W['ctrl'][0], W['lam_r']
        x2, y2 = self._xy(B, W)
        Z = O[1]['target']
        lam = mask(B, ctrl, ite(B, f1, lam_r, Z))
        if self.patched:
            f1o = B.XOR(f1, B.AND(ctrl, _s(is_zero(B, x2))))
        else:
            anc = eq(B, Z, lam_r)
            f1o = B.XOR(f1, B.AND(anc ^ 1, eq(B, lam, lam_r)))
        return dict(f1=[f1o], ctrl=[ctrl], a=W['a'], b=W['b'], x=x2, y=y2, lam=lam, lam_r=lam_r)


class Step3(_Base):
    name = 'step3'

    def pre(self, B, W, O):
        return self.rng(B, W['a'], W['b'], W['x'], W['y'], W['lam'])

    def calls(self, B, W, O):
        return [('MM', {'x': W['x'], 'y': W['lam']})]

    def splits(self, B, W, O):
        """Исходы трёх приведений по модулю: 2a>=p; a+(2a mod p)>=p; x+(3a mod p)>=p; и y<Z1 (CModSub)."""
        p, a, x = self.p, W['a'], W['x']
        P = lambda w: const(p, w)
        k1 = lambda B, W, O: lt(B, add(B, a, a), P(self.n + 1)) ^ 1
        k2 = lambda B, W, O: lt(B, add(B, a, self.N['T1']), P(self.n + 1)) ^ 1          # a + T1 >= p
        k3 = lambda B, W, O: lt(B, add(B, x, self.N['T2']), P(self.n + 1)) ^ 1          # x + T2 >= p
        k4 = lambda B, W, O: lt(B, W['y'], O[0]['target'])
        out = [('2a>=p', k1), ('a+2a>=p', k2)]
        if W['ctrl'][0] != ZERO:        # при ctrl=0 управляемые сложения тождественны: делить незачем
            out += [('x+3a>=p', k3), ('y<Z1', k4)]
        return out

    # Именованные промежуточные слова (определения — в контрольных точках): T1 = 2a mod p, T2 = (a + T1) mod p = 3a mod p.
    named = ('T1', 'T2')

    def definitions(self, B, W):
        """Определения именованных слов (те же формулы, что в контрольных точках)."""
        T1 = moddbl(B, W['a'], self.p)
        return dict(T1=T1, T2=modadd(B, W['a'], T1, self.p))

    def checkpoints(self, prog):
        from .core import region_end, region_outputs
        p = self.p
        out = lambda t, k, r: region_outputs(prog, t, k, r)
        return [(region_end(prog, 'ModDbl', 0), out('ModDbl', 0, 'x'), lambda B, W, O, N: moddbl(B, W['a'], p), 'T1'),
                (region_end(prog, 'ModAdd', 0), out('ModAdd', 0, 'y'), lambda B, W, O, N: modadd(B, W['a'], N['T1'], p), 'T2'),
                (region_end(prog, 'Adjoint', 0), out('Adjoint', 0, 'y'), lambda B, W, O, N: N['T1'], None),   # ModAdd^-1
                (region_end(prog, 'Adjoint', 1), out('Adjoint', 1, 'x'), lambda B, W, O, N: W['a'], None)]    # ModDbl^-1

    def post(self, B, W, O):
        p, ctrl, a = self.p, W['ctrl'][0], W['a']
        Z1 = O[0]['target']
        y3 = ite(B, ctrl, modsub(B, W['y'], Z1, p), W['y'])
        # x3 = (x + 3a mod p) mod p, где 3a mod p = T2 (определение T2 доказано в контрольной точке)
        x3 = ite(B, ctrl, modadd(B, W['x'], self.N['T2'], p), W['x'])
        return dict(ctrl=[ctrl], a=a, b=W['b'], x=x3, y=y3, lam=W['lam'])


class Step4(_Base):
    name = 'step4'
    range_words = ('x', 'lam')

    def pre(self, B, W, O):
        return self.rng(B, W['x'], W['lam'])

    def calls(self, B, W, O):
        x4 = modsub(B, W['x'], O[0]['target'], self.p)
        return [('MM', {'x': W['lam'], 'y': W['lam']}), ('MM', {'x': x4, 'y': W['lam']})]

    def post(self, B, W, O):
        x4 = modsub(B, W['x'], O[0]['target'], self.p)
        return dict(x=x4, y=xor(B, W['y'], O[1]['target']), lam=W['lam'])


class Step5(_Base):
    """Предусловие P5 (освобождение lam): lam = ctrl*(Z xor [K=0]*lam_r), K = KINV(x), Z = MM(K, y).
    Его выполнение доказывается на уровне ECAdd (E2)."""
    name = 'step5'
    range_words = ('a', 'b', 'x', 'y', 'lam_r')

    def pre(self, B, W, O):
        ctrl, lam_r = W['ctrl'][0], W['lam_r']
        K, Z = O[0]['x'], O[1]['target']
        p5 = eq(B, W['lam'], mask(B, ctrl, xor(B, Z, mask(B, is_zero(B, K), lam_r))))
        r = self.rng(B, W['a'], W['b'], W['x'], W['y'], W['lam_r'])
        return r if self.observe_on else r + [p5]

    def observe(self, B, W, O):
        """Вариант без P5: конечное значение проводов lam (их освобождает Free) = буквальная формула."""
        ctrl, lam_r, K, Z = W['ctrl'][0], W['lam_r'], O[0]['x'], O[1]['target']
        return dict(lam=xor(B, W['lam'], mask(B, ctrl, xor(B, Z, mask(B, is_zero(B, K), lam_r)))))

    def calls(self, B, W, O):
        return [('KINV', {'x': W['x']}), ('MM', {'x': O[0]['x'], 'y': W['y']})]

    def splits(self, B, W, O):
        p, ctrl = self.p, W['ctrl'][0]
        xn = lambda: ite(B, ctrl, modneg(B, W['x'], p), W['x'])
        out = [('x=0', is_zero(B, W['x'])),
               ('a+xn>=p', lambda B, W, O: lt(B, add(B, W['a'], xn()), const(p, self.n + 1)) ^ 1)]
        if ctrl != ZERO:
            out.append(('y<b', lt(B, W['y'], W['b'])))
        return out

    def post(self, B, W, O):
        p, ctrl = self.p, W['ctrl'][0]
        xn = ite(B, ctrl, modneg(B, W['x'], p), W['x'])
        x5 = modadd(B, W['a'], xn, p)
        y5 = ite(B, ctrl, modsub(B, W['y'], W['b'], p), W['y'])
        return dict(ctrl=[ctrl], a=W['a'], b=W['b'], x=x5, y=y5, lam_r=W['lam_r'])


class Step6(_Base):
    range_words = ('a', 'b')
    """Буквальное предусловие P6: все пять флагов в конце = 0 и аргументы CModSub/CModAdd < p.
    Его выполнение доказывается на уровне ECAdd (E2)."""
    name = 'step6'

    def literal(self, B, W):
        p = self.p
        f1, f2, f3, f4, ctrl = (W[k][0] for k in ('f1', 'f2', 'f3', 'f4', 'ctrl'))
        a, b, x, y = W['a'], W['b'], W['x'], W['y']
        za, zb, zx, zy = (_s(is_zero(B, w)) for w in (a, b, x, y))
        if self.patched:
            nn = B.AND(f3 ^ 1, f4 ^ 1)
            ctrl1 = B.XOR(B.XOR(ctrl, nn), B.AND(B.AND(f1, f2), nn))
        else:
            ctrl1 = B.XOR(ctrl, all_of(B, f2 ^ 1, f3 ^ 1, f4 ^ 1))
        f1a = _s(B.XOR(B.XOR(f1, all_of(B, za, zx, zy)), all_of(B, zx, za, zb)))
        f2a = _s(B.XOR(B.XOR(f2, all_of(B, za, zb, zy)), all_of(B, zx, zy, zb)))
        x1 = xor(B, x, mask(B, f4, a))
        y1 = xor(B, y, mask(B, f4, b))
        anc = _s(B.AND(f1a, f2a))
        x2 = ite(B, anc, modsub(B, x1, a, p), x1)
        y2 = ite(B, anc, modadd(B, y1, b, p), y1)
        f4a = _s(B.XOR(f4, B.AND(eq(B, a, x2), eq(B, b, y2))))
        f3a = _s(B.XOR(f3, B.AND(za, zb)))
        m = _s(B.AND(is_zero(B, x2), is_zero(B, y2)))
        f1b, f2b = B.XOR(f1a, m), B.XOR(f2a, m)
        return dict(x1=x1, y1=y1, x2=x2, y2=y2, anc=anc, flags=(ctrl1, f1b, f2b, f3a, f4a))

    def flag_prune(self, B, W):
        """Условие ctrl=0 в конце зависит только от флагов: если при данных флагах оно ложно — случай пуст."""
        f1, f2, f3, f4, ctrl = (W[k][0] for k in ('f1', 'f2', 'f3', 'f4', 'ctrl'))
        if self.patched:
            nn = B.AND(f3 ^ 1, f4 ^ 1)
            ctrl1 = B.XOR(B.XOR(ctrl, nn), B.AND(B.AND(f1, f2), nn))
        else:
            ctrl1 = B.XOR(ctrl, all_of(B, f2 ^ 1, f3 ^ 1, f4 ^ 1))
        return ctrl1 ^ 1 if not self.observe_on else ONE

    def splits(self, B, W, O):
        """Нули слов (управляют флагами и anc) и исходы приведений в CModSub/CModAdd."""
        p = self.p
        out = [(f'{k}=0', is_zero(B, W[k]), k) for k in ('a', 'b', 'x', 'y')]   # при «=1» слово заменяется нулями
        def k5(B, W, O):
            L = self.literal(B, W)
            return ZERO if L['anc'] == ZERO else lt(B, L['x1'], W['a'])          # без CModSub делить незачем
        def k6(B, W, O):
            L = self.literal(B, W)
            return ZERO if L['anc'] == ZERO else lt(B, add(B, L['y1'], W['b']), const(p, self.n + 1)) ^ 1
        out.append(('x1<a', k5))
        out.append(('y1+b>=p', k6))
        return out

    def pre(self, B, W, O):
        L = self.literal(B, W)
        r = self.rng(B, W['a'], W['b'], L['x1'], L['y1'])
        return r if self.observe_on else r + [f ^ 1 for f in L['flags']]

    def observe(self, B, W, O):
        """Вариант без условия на флаги: конечные значения пяти флагов (их освобождает Free) = буквальные формулы."""
        ctrl1, f1b, f2b, f3a, f4a = self.literal(B, W)['flags']
        return dict(ctrl=[ctrl1], f1=[f1b], f2=[f2b], f3=[f3a], f4=[f4a])

    def post(self, B, W, O):
        L = self.literal(B, W)
        return dict(a=W['a'], b=W['b'], x=L['x2'], y=L['y2'])


STEPS = {1: Step1, 2: Step2, 3: Step3, 4: Step4, 5: Step5, 6: Step6}
