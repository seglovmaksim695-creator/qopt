# SPDX-License-Identifier: Apache-2.0
"""Registry of contracts for Qualtran bloqs: primitives (proved directly) and composite specifications.

``contract_for(bloq)`` returns a ``Contract`` or ``None``.  Kind ``primitive``: small arithmetic building blocks
(addition, subtraction, comparison, copy, XOR, swap, table lookup, basis states ...) proved on their own operations
by kprove.  Kind ``spec``: composite blocks proved from their children (``specs``).  Every entry is checked against
the classical simulation of Qualtran (``bloq.call_classically``) on all inputs of small instances by ``selfcheck``.

Only unsigned types (``QUInt``, ``QMontgomeryUInt``, ``BQUInt``, ``QAny`` where the operation is bitwise) get a
contract; signed or symbolic parameters give ``None`` (the node is then derived from its children, or UNKNOWN).
"""
from __future__ import annotations

from . import expr as X
from . import specs as S
from .contract import Contract

__all__ = ['contract_for', 'selfcheck', 'SELFCHECK_INSTANCES', 'to_classical', 'from_classical']

V = X.var


def _int(v):
    try:
        iv = int(v)
    except Exception:
        return None
    return iv if iv == v else None


def _unsigned(dt):
    from qualtran import QUInt, QMontgomeryUInt, BQUInt, QAny, QBit
    return type(dt) in (QUInt, QMontgomeryUInt, BQUInt, QAny, QBit)


def _prim(name, inputs, outputs, out, pre=1, defs=(), note=''):
    return Contract(inputs=inputs, outputs=outputs, out=out, pre=pre, defs=list(defs), name=name, kind='primitive',
                    note=note).validate()


def _bits(dt):
    n = _int(dt.num_qubits)
    return n if n and n > 0 else None


def contract_for(bloq):
    """the registry contract of ``bloq`` or None"""
    try:
        return _contract_for(bloq)
    except S.SpecDomainError:
        return None


def _contract_for(bloq):
    from qualtran.bloqs.basic_gates import (XGate, CNOT, Toffoli, Swap, CSwap, TwoBitSwap, TwoBitCSwap, ZeroState,
                                            OneState, ZeroEffect, OneEffect, IntState, IntEffect)
    from qualtran.bloqs.mcmt import And, MultiTargetCNOT
    from qualtran.bloqs.arithmetic import (Add, AddK, LessThanConstant, Xor, XorK, LessThanEqual, GreaterThan,
                                           Equals, EqualsAConstant, BitwiseNot)
    from qualtran.bloqs.arithmetic.comparison import LinearDepthGreaterThan, CLinearDepthGreaterThan
    from qualtran.bloqs.arithmetic.controlled_addition import CAdd
    from qualtran.bloqs.data_loading.qroam_clean import QROAMClean, QROAMCleanAdjoint, QROAMCleanAdjointWrapper
    from qualtran.bloqs.mod_arithmetic import ModAdd, CModAdd, ModSub, CModSub, ModNeg, CModNeg, ModDbl
    from qualtran.bloqs.mod_arithmetic.mod_multiplication import DirtyOutOfPlaceMontgomeryModMul
    from qualtran import Controlled
    t = type(bloq)
    # ---------------------------------------------------------------- single gates and basis states
    if t is XGate:
        return _prim('X', {'q': 1}, {'q': 1}, {'q': X.sub(1, V('q'))})
    if t is CNOT:
        return _prim('CNOT', {'ctrl': 1, 'target': 1}, {'ctrl': 1, 'target': 1},
                     {'ctrl': V('ctrl'), 'target': X.bxor(V('target'), V('ctrl'))})
    if t is Toffoli:
        return _prim('Toffoli', {'ctrl': 2, 'target': 1}, {'ctrl': 2, 'target': 1},
                     {'ctrl': V('ctrl'), 'target': X.bxor(V('target'), X.eq(V('ctrl'), 3))})
    if t is And:
        c1, c2 = _int(bloq.cv1), _int(bloq.cv2)
        if c1 not in (0, 1) or c2 not in (0, 1):
            return None
        act = X.eq(V('ctrl'), 2 * c1 + c2)
        if bloq.uncompute:
            return _prim(f'And^-1({c1},{c2})', {'ctrl': 2, 'target': 1}, {'ctrl': 2}, {'ctrl': V('ctrl')},
                         pre=X.eq(V('target'), act))
        return _prim(f'And({c1},{c2})', {'ctrl': 2}, {'ctrl': 2, 'target': 1}, {'ctrl': V('ctrl'), 'target': act})
    if t in (ZeroState, OneState):
        return _prim(t.__name__, {}, {'q': 1}, {'q': 1 if t is OneState else 0})
    if t in (ZeroEffect, OneEffect):
        return _prim(t.__name__, {'q': 1}, {}, {}, pre=X.eq(V('q'), 1 if t is OneEffect else 0))
    if t is IntState:
        n, k = _int(bloq.bitsize), _int(bloq.val)
        if not n or k is None or not 0 <= k < (1 << n):
            return None
        return Contract({}, {'val': n}, {'val': k}, name=f'IntState({k}, {n})', kind='spec').validate()
    if t is IntEffect:
        n, k = _int(bloq.bitsize), _int(bloq.val)
        if not n or k is None or not 0 <= k < (1 << n):
            return None
        return Contract({'val': n}, {}, {}, pre=X.eq(V('val'), k), name=f'IntEffect({k}, {n})', kind='spec').validate()
    if t in (TwoBitSwap,):
        return _prim('TwoBitSwap', {'x': 1, 'y': 1}, {'x': 1, 'y': 1}, {'x': V('y'), 'y': V('x')})
    if t in (TwoBitCSwap,):
        c = X.eq(V('ctrl'), 1)
        return _prim('TwoBitCSwap', {'ctrl': 1, 'x': 1, 'y': 1}, {'ctrl': 1, 'x': 1, 'y': 1},
                     {'ctrl': V('ctrl'), 'x': X.ite(c, V('y'), V('x')), 'y': X.ite(c, V('x'), V('y'))})
    if t is Swap:
        n = _int(bloq.bitsize)
        return _prim(f'Swap({n})', {'x': n, 'y': n}, {'x': n, 'y': n}, {'x': V('y'), 'y': V('x')}) if n else None
    if t is CSwap:
        n = _int(bloq.bitsize)
        if not n:
            return None
        c = X.eq(V('ctrl'), 1)
        return _prim(f'CSwap({n})', {'ctrl': 1, 'x': n, 'y': n}, {'ctrl': 1, 'x': n, 'y': n},
                     {'ctrl': V('ctrl'), 'x': X.ite(c, V('y'), V('x')), 'y': X.ite(c, V('x'), V('y'))})
    if t is MultiTargetCNOT:
        n = _int(bloq.bitsize)
        if not n:
            return None
        return _prim(f'MultiTargetCNOT({n})', {'control': 1, 'targets': n}, {'control': 1, 'targets': n},
                     {'control': V('control'), 'targets': X.bxor(V('targets'), X.ite(X.eq(V('control'), 1), (1 << n) - 1, 0))})
    # ---------------------------------------------------------------- arithmetic primitives
    if t is Add:
        if not (_unsigned(bloq.a_dtype) and _unsigned(bloq.b_dtype)):
            return None
        na, nb = _bits(bloq.a_dtype), _bits(bloq.b_dtype)
        if not na or not nb:
            return None
        return _prim(f'Add({na},{nb})', {'a': na, 'b': nb}, {'a': na, 'b': nb},
                     {'a': V('a'), 'b': X.mod(X.add(V('a'), V('b')), 1 << nb)})
    if t is AddK:
        n, k = _bits(bloq.dtype), _int(bloq.k)
        if not _unsigned(bloq.dtype) or not n or k is None:
            return None
        return _prim(f'AddK({n},{k})', {'x': n}, {'x': n}, {'x': X.mod(X.add(V('x'), k), 1 << n)})
    if t is CAdd:
        if not (_unsigned(bloq.a_dtype) and _unsigned(bloq.b_dtype)):
            return None
        na, nb, cv = _bits(bloq.a_dtype), _bits(bloq.b_dtype), _int(bloq.cv)
        if not na or not nb or cv not in (0, 1):
            return None
        return _prim(f'CAdd({na},{nb},cv={cv})', {'ctrl': 1, 'a': na, 'b': nb}, {'ctrl': 1, 'a': na, 'b': nb},
                     {'ctrl': V('ctrl'), 'a': V('a'),
                      'b': X.ite(X.eq(V('ctrl'), cv), X.mod(X.add(V('a'), V('b')), 1 << nb), V('b'))})
    if t is LessThanConstant:
        n, k = _int(bloq.bitsize), _int(bloq.less_than_val)
        if not n or k is None:
            return None
        return _prim(f'LessThanConstant({n},{k})', {'x': n, 'target': 1}, {'x': n, 'target': 1},
                     {'x': V('x'), 'target': X.bxor(V('target'), X.lt(V('x'), k))})
    if t is EqualsAConstant:
        n, k = _int(bloq.bitsize), _int(bloq.val)
        if not n or k is None:
            return None
        return _prim(f'EqualsAConstant({n},{k})', {'x': n, 'target': 1}, {'x': n, 'target': 1},
                     {'x': V('x'), 'target': X.bxor(V('target'), X.eq(V('x'), k))})
    if t is LessThanEqual:
        nx, ny = _int(bloq.x_bitsize), _int(bloq.y_bitsize)
        if not nx or not ny:
            return None
        return _prim(f'LessThanEqual({nx},{ny})', {'x': nx, 'y': ny, 'target': 1}, {'x': nx, 'y': ny, 'target': 1},
                     {'x': V('x'), 'y': V('y'), 'target': X.bxor(V('target'), X.le(V('x'), V('y')))})
    if t is GreaterThan:
        na, nb = _int(bloq.a_bitsize), _int(bloq.b_bitsize)
        if not na or not nb:
            return None
        return _prim(f'GreaterThan({na},{nb})', {'a': na, 'b': nb, 'target': 1}, {'a': na, 'b': nb, 'target': 1},
                     {'a': V('a'), 'b': V('b'), 'target': X.bxor(V('target'), X.lt(V('b'), V('a')))})
    if t is LinearDepthGreaterThan:
        n = _int(bloq.bitsize)
        if not n or bloq.signed:
            return None
        return _prim(f'LinearDepthGreaterThan({n})', {'a': n, 'b': n, 'target': 1}, {'a': n, 'b': n, 'target': 1},
                     {'a': V('a'), 'b': V('b'), 'target': X.bxor(V('target'), X.lt(V('b'), V('a')))})
    if t is CLinearDepthGreaterThan:
        n, cv = _bits(bloq.dtype), _int(bloq.cv)
        if not n or not _unsigned(bloq.dtype) or cv not in (0, 1):
            return None
        return _prim(f'CLinearDepthGreaterThan({n},cv={cv})', {'ctrl': 1, 'a': n, 'b': n, 'target': 1},
                     {'ctrl': 1, 'a': n, 'b': n, 'target': 1},
                     {'ctrl': V('ctrl'), 'a': V('a'), 'b': V('b'),
                      'target': X.bxor(V('target'), X.and_(X.eq(V('ctrl'), cv), X.lt(V('b'), V('a'))))})
    if t is Equals:
        n = _bits(bloq.dtype)
        if not n:
            return None
        return _prim(f'Equals({n})', {'x': n, 'y': n, 'target': 1}, {'x': n, 'y': n, 'target': 1},
                     {'x': V('x'), 'y': V('y'), 'target': X.bxor(V('target'), X.eq(V('x'), V('y')))})
    if t is Xor:
        n = _bits(bloq.dtype)
        if not n:
            return None
        return _prim(f'Xor({n})', {'x': n, 'y': n}, {'x': n, 'y': n}, {'x': V('x'), 'y': X.bxor(V('x'), V('y'))})
    if t is XorK:
        n, k = _bits(bloq.dtype), _int(bloq.k)
        if not n or k is None or not 0 <= k < (1 << n):
            return None
        return _prim(f'XorK({n},{k})', {'x': n}, {'x': n}, {'x': X.bxor(V('x'), k)})
    if t is BitwiseNot:
        n = _bits(bloq.dtype)
        if not n or not _unsigned(bloq.dtype):
            return None
        return _prim(f'BitwiseNot({n})', {'x': n}, {'x': n}, {'x': X.sub((1 << n) - 1, V('x'))})
    # ---------------------------------------------------------------- table lookup (rule C1 for the adjoint)
    if t in (QROAMClean, QROAMCleanAdjoint, QROAMCleanAdjointWrapper):
        return _qroam(bloq)
    # ---------------------------------------------------------------- controlled primitives, adjoints, Always
    if t is Controlled:
        return _controlled(bloq)
    from qualtran import Adjoint
    if t is Adjoint:
        return _adjoint(bloq.subbloq)
    if t.__name__ == 'Always' and hasattr(bloq, 'subbloq'):
        c = _contract_for(bloq.subbloq)
        return c if c is not None and c.kind == 'primitive' else None
    # ---------------------------------------------------------------- composite specifications
    if t in (ModAdd,):
        n, p = _int(bloq.bitsize), _int(bloq.mod)
        return S.mod_add_spec(n, p) if n and p else None
    if t is CModAdd:
        n, p, cv = _bits(bloq.dtype), _int(bloq.mod), _int(bloq.cv)
        return S.cmod_add_spec(n, p, cv) if n and p and _unsigned(bloq.dtype) else None
    if t is ModSub:
        n, p = _bits(bloq.dtype), _int(bloq.mod)
        return S.mod_sub_spec(n, p) if n and p and _unsigned(bloq.dtype) else None
    if t is CModSub:
        n, p, cv = _bits(bloq.dtype), _int(bloq.mod), _int(getattr(bloq, 'cv', 1))
        return S.cmod_sub_spec(n, p, cv) if n and p and _unsigned(bloq.dtype) else None
    if t is ModNeg:
        n, p = _bits(bloq.dtype), _int(bloq.mod)
        return S.mod_neg_spec(n, p) if n and p and _unsigned(bloq.dtype) else None
    if t is CModNeg:
        n, p, cv = _bits(bloq.dtype), _int(bloq.mod), _int(getattr(bloq, 'cv', 1))
        return S.cmod_neg_spec(n, p, cv) if n and p and _unsigned(bloq.dtype) else None
    if t is ModDbl:
        n, p = _bits(bloq.dtype), _int(bloq.mod)
        return S.mod_dbl_spec(n, p) if n and p and _unsigned(bloq.dtype) else None
    if t is DirtyOutOfPlaceMontgomeryModMul:
        if bloq.uncompute:
            return None
        n, w, p = _int(bloq.bitsize), _int(bloq.window_size), _int(bloq.mod)
        return S.montgomery_spec(n, w, p) if n and w and p else None
    return None


def _qroam(bloq):
    from qualtran.bloqs.data_loading.qroam_clean import QROAMClean
    try:
        if not bloq.has_data():
            return None
        sel = list(bloq.selection_registers)
        tgt = list(bloq.target_registers)
        ctl = list(getattr(bloq, 'control_registers', ()))
        datas = bloq.data_or_shape
    except Exception:
        return None
    if len(sel) != 1 or len(tgt) != 1 or ctl or len(datas) != 1:
        return None
    s, tr = sel[0], tgt[0]
    if s.shape or tuple(tr.shape) not in ((), (1,)):
        return None
    tname = tr.name if not tr.shape else tr.name + '[0]'
    import numpy as np
    d = np.asarray(datas[0], dtype=object).reshape(-1)
    table = [_int(v) for v in d]
    nt, ns = _int(tr.dtype.num_qubits), _int(s.dtype.num_qubits)
    if any(v is None or v < 0 or v >> nt for v in table) or len(table) > (1 << ns):
        return None
    if len(table) < (1 << ns):
        return None       # a selection outside the table: no contract (the iteration length would be a precondition)
    names = {r.name for r in bloq.signature}
    if type(bloq) is QROAMClean:
        if names != {s.name, tr.name}:
            return None   # junk registers (block size > 1): no contract
        return _prim(f'QROAMClean(table of {len(table)})', {s.name: ns}, {s.name: ns, tname: nt},
                     {s.name: V(s.name), tname: X.lut(V(s.name), table)})
    if names != {s.name, tr.name}:
        return None
    return _prim(f'QROAMCleanAdjoint(table of {len(table)})', {s.name: ns, tname: nt}, {s.name: ns},
                 {s.name: V(s.name)}, pre=X.eq(V(tname), X.lut(V(s.name), table)),
                 note='measurement-based uncomputation (bridge rule C1, assumptions A1-A4)')


_SELF_INVERSE = ('XGate', 'CNOT', 'Toffoli', 'Swap', 'CSwap', 'TwoBitSwap', 'TwoBitCSwap', 'MultiTargetCNOT', 'Xor',
                 'XorK', 'BitwiseNot', 'LessThanConstant', 'LessThanEqual', 'GreaterThan', 'LinearDepthGreaterThan',
                 'CLinearDepthGreaterThan', 'Equals', 'EqualsAConstant')


def _adjoint(sub):
    """the inverse of a primitive: self-inverse gates and comparisons (an XOR into the target) keep their contract;
    additions get the subtraction; Controlled and Always are unwrapped.  Checked against call_classically."""
    from qualtran import Controlled
    from qualtran.bloqs.arithmetic import Add, AddK
    from qualtran.bloqs.arithmetic.controlled_addition import CAdd
    t = type(sub)
    if t.__name__ == 'Always' and hasattr(sub, 'subbloq'):
        return _adjoint(sub.subbloq)
    if t.__name__ in _SELF_INVERSE:
        c = _contract_for(sub)
        if c is None or c.inputs != c.outputs:
            return None
        return _prim(f'Adjoint[{c.name}]', c.inputs, c.outputs, c.out, pre=c.pre, defs=c.defs)
    if t is AddK:
        n, k = _bits(sub.dtype), _int(sub.k)
        if not _unsigned(sub.dtype) or not n or k is None:
            return None
        return _prim(f'Adjoint[AddK({n},{k})]', {'x': n}, {'x': n}, {'x': X.mod(X.sub(V('x'), k), 1 << n)})
    if t is Add:
        if not (_unsigned(sub.a_dtype) and _unsigned(sub.b_dtype)):
            return None
        na, nb = _bits(sub.a_dtype), _bits(sub.b_dtype)
        if not na or not nb:
            return None
        return _prim(f'Adjoint[Add({na},{nb})]', {'a': na, 'b': nb}, {'a': na, 'b': nb},
                     {'a': V('a'), 'b': X.mod(X.sub(V('b'), V('a')), 1 << nb)})
    if t is CAdd:
        if not (_unsigned(sub.a_dtype) and _unsigned(sub.b_dtype)):
            return None
        na, nb, cv = _bits(sub.a_dtype), _bits(sub.b_dtype), _int(sub.cv)
        if not na or not nb or cv not in (0, 1):
            return None
        return _prim(f'Adjoint[CAdd({na},{nb},cv={cv})]', {'ctrl': 1, 'a': na, 'b': nb}, {'ctrl': 1, 'a': na, 'b': nb},
                     {'ctrl': V('ctrl'), 'a': V('a'),
                      'b': X.ite(X.eq(V('ctrl'), cv), X.mod(X.sub(V('b'), V('a')), 1 << nb), V('b'))})
    if t is Controlled:
        inner = _adjoint(sub.subbloq)
        if inner is None:
            return None
        return _controlled(sub, sub_contract=inner, prefix='Adjoint')
    return None


def _controlled(bloq, sub_contract=None, prefix=''):
    """Controlled(sub) with one QBit control: the sub-contract when the control is active, identity otherwise.
    Only for a sub-contract whose outputs are exactly its inputs (THRU registers of equal width)."""
    from qualtran import QBit
    spec = bloq.ctrl_spec
    try:
        qdts = list(spec.qdtypes)
        cvs = [_int(c) for c in spec.cvs]
        shapes = [tuple(getattr(c, 'shape', ())) for c in spec.cvs]
    except Exception:
        return None
    if len(qdts) != 1 or type(qdts[0]) is not QBit or cvs[0] not in (0, 1) or shapes[0] != ():
        return None
    sub = sub_contract if sub_contract is not None else _contract_for(bloq.subbloq)
    if sub is None or sub.kind != 'primitive' or sub.inputs != sub.outputs or sub.folds:
        return None
    names = [r.name for r in bloq.signature]
    ctrl_names = [n for n in names if n not in sub.inputs]
    if len(ctrl_names) != 1:
        return None
    c = ctrl_names[0]
    if c in sub.inputs:
        return None
    act = X.eq(V(c), cvs[0])
    out = {c: V(c)}
    for r in sub.outputs:
        out[r] = X.ite(act, sub.out[r], V(r))
    pre = sub.pre if sub.pre == 1 else X.or_(X.not_(act), sub.pre)
    inputs = {c: 1, **sub.inputs}
    return _prim(f'{prefix}Controlled[{sub.name}, cv={cvs[0]}]', inputs, dict(inputs), out, pre=pre, defs=sub.defs)


# ------------------------------------------------------------------------------------------------ self-check
def to_classical(reg, packed):
    import numpy as np
    shape = tuple(int(s) for s in reg.shape)
    if not shape:
        return packed
    if len(shape) == 1 and int(reg.bitsize) == 1:
        m = shape[0]
        return np.array([(packed >> (m - 1 - k)) & 1 for k in range(m)], dtype=object)
    raise ValueError('element-wise registers are not packed')


def from_classical(reg, value):
    shape = tuple(int(s) for s in reg.shape)
    if not shape:
        return int(value)
    m = shape[0]
    return sum(int(value[k]) << (m - 1 - k) for k in range(m))


def _small_instances():
    """small instances of every registry type: (bloq, domain(values) -> bool or None, kind of Qualtran or None)"""
    import numpy as np
    from qualtran import QUInt, QMontgomeryUInt
    from qualtran.bloqs.basic_gates import (XGate, CNOT, Toffoli, Swap, CSwap, TwoBitSwap, TwoBitCSwap, IntState,
                                            IntEffect)
    from qualtran.bloqs.mcmt import And, MultiTargetCNOT
    from qualtran.bloqs.arithmetic import (Add, AddK, LessThanConstant, Xor, XorK, LessThanEqual, GreaterThan,
                                           Equals, EqualsAConstant, BitwiseNot)
    from qualtran.bloqs.arithmetic.comparison import LinearDepthGreaterThan, CLinearDepthGreaterThan
    from qualtran.bloqs.arithmetic.controlled_addition import CAdd
    from qualtran.bloqs.data_loading.qroam_clean import QROAMClean
    from qualtran.bloqs.mod_arithmetic import ModAdd, CModAdd, ModSub, CModSub, ModNeg, CModNeg, ModDbl
    out = [XGate(), CNOT(), Toffoli(), TwoBitSwap(), TwoBitCSwap(), Swap(2), CSwap(2), MultiTargetCNOT(3)]
    out += [And(a, b) for a in (0, 1) for b in (0, 1)] + [And(1, 0, uncompute=True), And(1, 1, uncompute=True)]
    out += [IntState(5, 3), IntEffect(5, 3), IntState(0, 2)]
    out += [Add(QUInt(3)), Add(QUInt(2), QUInt(3)), Add(QMontgomeryUInt(3)), AddK(QUInt(3), 5), AddK(QUInt(3), -5),
            AddK(QUInt(4), -11), AddK(QUInt(3), 2).controlled(), AddK(QUInt(3), -3).controlled(),
            CAdd(QUInt(2), QUInt(3)), CAdd(QUInt(3), QUInt(3), cv=0), CAdd(QMontgomeryUInt(2), QMontgomeryUInt(3))]
    out += [LessThanConstant(3, 5), LessThanConstant(3, 0), EqualsAConstant(3, 5),
            LessThanEqual(3, 3), GreaterThan(3, 3), LinearDepthGreaterThan(3), CLinearDepthGreaterThan(QUInt(3)),
            CLinearDepthGreaterThan(QUInt(3), cv=0), Equals(QUInt(3)), Xor(QUInt(3)), XorK(QUInt(3), 5),
            BitwiseNot(QUInt(3))]
    out += [AddK(QUInt(3), 5).adjoint(), AddK(QUInt(4), -11).adjoint(), Add(QUInt(3)).adjoint(),
            Add(QUInt(2), QUInt(3)).adjoint(), CAdd(QUInt(2), QUInt(3)).adjoint(), CAdd(QUInt(3), QUInt(3), cv=0).adjoint(),
            Xor(QUInt(3)).adjoint(), XorK(QUInt(3), 5).adjoint(), BitwiseNot(QUInt(3)).adjoint(),
            LessThanConstant(3, 5).adjoint(), AddK(QUInt(3), 2).controlled().adjoint(), CNOT().adjoint()]
    out += [QROAMClean([np.array([3, 10], dtype=object)], selection_bitsizes=(1,), target_bitsizes=(4,)),
            QROAMClean([np.array([0, 9, 2, 15], dtype=object)], selection_bitsizes=(2,), target_bitsizes=(4,))]
    out += [ModAdd(3, 7), ModAdd(4, 13), CModAdd(QUInt(3), 7), CModAdd(QUInt(4), 11, cv=0), ModSub(QUInt(3), 7),
            CModSub(QUInt(3), 5), ModNeg(QUInt(3), 7), CModNeg(QUInt(3), 5), ModDbl(QUInt(4), 13)]
    return out


def selfcheck(instances=None, *, max_inputs=1 << 14):
    """For every instance with a registry contract: enumerate every input of its domain and compare the contract with
    ``bloq.call_classically``.  Returns a list of dicts (``ok``, ``checked``, ``first_mismatch``)."""
    import itertools
    rows = []
    for bloq in (instances if instances is not None else _small_instances()):
        c = contract_for(bloq)
        row = dict(bloq=repr(bloq)[:120], contract=c.name if c else None, ok=False, checked=0, first_mismatch=None)
        if c is None:
            row['why'] = 'no contract'
            rows.append(row)
            continue
        lefts = list(bloq.signature.lefts())
        rights = list(bloq.signature.rights())
        names = [r.name for r in lefts]
        doms = [range(1 << c.inputs[n]) for n in names]
        total = 1
        for d in doms:
            total *= len(d)
        if total > max_inputs:
            row['why'] = f'{total} inputs > {max_inputs}'
            rows.append(row)
            continue
        bad = None
        cnt = 0
        for vals in itertools.product(*doms):
            v = dict(zip(names, vals))
            if not c.in_domain(v):
                continue
            want = c.evaluate(v)
            try:
                res = bloq.call_classically(**{r.name: to_classical(r, v[r.name]) for r in lefts})
            except Exception as e:  # noqa: BLE001
                bad = dict(input=v, why=f'call_classically raised {type(e).__name__}: {e}')
                break
            if not isinstance(res, tuple):
                res = (res,)
            got = {r.name: from_classical(r, val) for r, val in zip(rights, res)}
            cnt += 1
            if got != want:
                bad = dict(input=v, contract=want, qualtran=got)
                break
        row.update(checked=cnt, ok=bad is None and cnt > 0, first_mismatch=bad)
        rows.append(row)
    return rows
