# SPDX-License-Identifier: Apache-2.0
"""The compositional prover: a node is proved from the contracts of its children, recursively.

``Engine(program).prove(node, contract)`` returns a result dict whose ``status`` is ``PROVED`` only when every check
below holds; anything else (a failed check, a refuted obligation, a Z3 time-out, an unsupported construct) is
``NOT_PROVED`` / ``UNKNOWN`` with the reason.

How a node is proved
--------------------
* a **leaf** (a node whose contract comes from the primitive registry, or a node without children) is proved
  directly on its own operations (``leaf.prove_leaf``: kprove ROBDD + exact AND/ERASE domains + independent
  simulation);
* a **composite** node is proved by walking its children in the order of the source and keeping a symbolic state
  ``wire -> 0 | 1 | bit j of value v``, where every value is an integer expression over the inputs of the node
  (``expr``).  For every child:

  1. *coverage*: the children own consecutive, non-overlapping ranges of operations that tile the range of the node
     exactly (an operation of the node that belongs to no child is refused, AUD-C-01);
  2. *freshness*: every wire the child assumes to be 0 (its operations and outputs minus its inputs) is the constant 0
     in the state and was never freed before (AUD-C-01);
  3. *wire flow*: the inputs of the child are read from its real input wires (from the source graph) in the state;
     its contract is applied to exactly these values; its outputs are written to its real output wires; every other
     wire it touched is 0 afterwards (its proof says so) (AUD-D-01);
  4. *domain*: the precondition of the child, on the values read, follows from the facts (Z3, linear integer
     arithmetic; unsat of the negation) -- otherwise the node is not proved (AUD7: domain before counting);
  5. allocations must be fresh, freed wires must provably be 0;
  then the outputs of the node must equal its contract and every other wire must be 0 (Z3).
* a node **without a contract** gets a *derived* contract (its outputs as expressions of its inputs; the
  preconditions of its children that do not follow from the facts of the node become its precondition); the
  parent must then discharge them.  The root always needs a contract.
* a **loop** (a run of children of the type ``Fold.body``) is accepted only with a ``Fold`` of the parent contract:
  every iteration is proved against ``Fold.body_contract``; its inputs must be the ones the fold states; the carried
  state after an iteration is a fresh symbol about which only the range and the invariant are known; the invariant is
  proved at iteration 0 and by one step lemma for a symbolic iteration.  Without a fold the children are unrolled.

Facts used by Z3: ranges of the inputs, the precondition of the node, the output ranges of every child (part of the
proof of that child), the invariants of fold symbols, and for derived contracts the deferred preconditions.
"""
from __future__ import annotations

import time

from . import expr as X
from .contract import Contract, ContractError
from .leaf import prove_leaf

__all__ = ['Engine', 'ProofFail']

ZERO_C = ('c', 0)
ONE_C = ('c', 1)


class ProofFail(Exception):
    def __init__(self, msg, status='NOT_PROVED', **detail):
        super().__init__(msg)
        self.status = status
        self.detail = detail


class _Val:
    __slots__ = ('expr', 'width')

    def __init__(self, expr, width):
        self.expr, self.width = expr, width


class _Comp:
    """the symbolic state of one composite node"""

    def __init__(self, engine, node, contract, z3):
        self.E, self.node, self.C, self.z3 = engine, node, contract, z3
        self.node_widths = node.widths_in()
        self.vals = []
        self.state = {}
        self.retired = set()
        self.touched = set()
        self.syms = {}
        self.n_obl = 0
        self.n_z3 = 0
        self.deferred = []
        self.derived = contract is None
        self.T = X.Z3Terms(z3, self._z3_leaf, self._leaf_bounds)
        self.var_terms = {}
        self.global_facts = []          # z3 Bools: precondition of the node, deferred conditions
        self.val_facts = {}             # vid -> [z3 Bool]   (the output range of the child that made it)
        self.sym_facts = {}             # sym name -> [z3 Bool]  (range, invariant of a fold)

    # ---------------------------------------------------------------- values and Z3 leaves
    def new_val(self, e, width):
        self.vals.append(_Val(e, width))
        vid = len(self.vals) - 1
        if width is not None:
            t = self.T.term(('val', vid))
            self.val_facts[vid] = [t >= 0, t < (1 << width)]
        return vid

    def add_sym(self, name, width):
        self.syms[name] = width
        t = self.T.term(('sym', name))
        self.sym_facts.setdefault(name, []).extend([t >= 0, t < (1 << width)])

    def add_fact(self, cond_expr):
        """a proved (or assumed-by-contract) 0/1 fact; attached to every symbol it mentions, else global"""
        c = self.T.cond(cond_expr)
        syms = X.refs_of(cond_expr, 'sym')
        if syms:
            for s in syms:
                self.sym_facts.setdefault(s, []).append(c)
        else:
            self.global_facts.append(c)

    def facts_for(self, e):
        """the facts about every value, symbol and input that e depends on (cone of influence)"""
        out = list(self.global_facts)
        seen_v, seen_s, seen_x = set(), set(), set()
        stack = [e]
        seen = set()
        while stack:
            x = stack.pop()
            if type(x) is int or id(x) in seen:
                continue
            seen.add(id(x))
            k = x[0]
            if k == 'val':
                vid = x[1]
                if vid not in seen_v:
                    seen_v.add(vid)
                    out.extend(self.val_facts.get(vid, ()))
                    stack.append(self.vals[vid].expr)
            elif k == 'sym':
                if x[1] not in seen_s:
                    seen_s.add(x[1])
                    out.extend(self.sym_facts.get(x[1], ()))
            elif k == 'var':
                if x[1] not in seen_x:
                    seen_x.add(x[1])
                    t = self.T.term(x)
                    out.extend([t >= 0, t < (1 << self.node_widths[x[1]])])
            elif k in ('ref', 'fold'):
                continue
            elif k == 'lut':
                stack.append(x[1])
            else:
                stack.extend(a for a in x[1:] if type(a) is tuple)
        # facts of symbols can mention further symbols (an invariant over several state registers)
        return out

    def _leaf_bounds(self, leaf):
        k = leaf[0]
        if k == 'var':
            return (0, (1 << self.node_widths[leaf[1]]) - 1)
        if k == 'val':
            v = self.vals[leaf[1]]
            if v.width is None:
                return self.T.bounds(v.expr)
            return (0, (1 << v.width) - 1)
        if k == 'sym':
            return (0, (1 << self.syms[leaf[1]]) - 1)
        raise X.ExprError(f'unexpected leaf {leaf!r}')

    def _z3_leaf(self, leaf):
        z3 = self.z3
        k = leaf[0]
        if k == 'var':
            t = self.var_terms.get(leaf[1])
            if t is None:
                t = self.var_terms[leaf[1]] = z3.Int('v_' + leaf[1])
            return t
        if k == 'val':
            return self.T.term(self.vals[leaf[1]].expr)
        if k == 'sym':
            return z3.Int('s_' + leaf[1])
        raise X.ExprError(f'unexpected leaf {leaf!r}')

    # ---------------------------------------------------------------- wires
    def read(self, wires):
        """the integer on these wires (LSB first) as an expression; runs of consecutive bits of one value become one
        slice (the whole value when the run covers it)"""
        parts, const, j, n = [], 0, 0, len(wires)
        st = self.state
        while j < n:
            s = st.get(wires[j], ZERO_C)
            if s[0] == 'c':
                const |= s[1] << j
                j += 1
                continue
            vid, k = s
            L = 1
            while j + L < n:
                t = st.get(wires[j + L], ZERO_C)
                if t[0] != vid or t[1] != k + L:
                    break
                L += 1
            w = self.vals[vid].width
            term = ('val', vid) if (k == 0 and w is not None and L >= w) else X.slc(('val', vid), k, L)
            parts.append(term if j == 0 else X.shl(term, j))
            j += L
        if const:
            parts.append(const)
        return X.sum_(parts) if parts else 0

    def write_value(self, wires, vid):
        for j, q in enumerate(wires):
            self.state[q] = (vid, j)

    # ---------------------------------------------------------------- obligations
    def obligation(self, cond, what, *, must=True):
        """cond: a 0/1 expression.  Proved valid under the facts, or (derived contracts, must=False) deferred."""
        self.n_obl += 1
        if cond == 1:
            return True
        if type(cond) is tuple and cond[0] == 'eq' and (cond[1] is cond[2] or cond[1] == cond[2]):
            return True
        z3 = self.z3
        try:
            if type(cond) is tuple and cond[0] == 'eq':
                ta, tb = self.T.term(cond[1]), self.T.term(cond[2])
                if ta.eq(tb) or z3.simplify(ta - tb).eq(z3.IntVal(0)):
                    return True
            c = self.T.cond(cond)
            if self.derived and not must:
                # a node without a contract: its open conditions go to its precondition (sound: the parent must
                # prove them); no solver call here (in context they are usually syntactic)
                self.deferred.append((cond, what))
                self.global_facts.append(c)
                return True
            facts = self.facts_for(cond)
        except X.ExprError as e:
            raise ProofFail(f'{what}: not expressible in linear integer arithmetic: {e}', status='UNKNOWN')
        s = z3.Solver()
        s.set('timeout', self.E.z3_timeout_ms)
        s.add(*facts)
        s.add(z3.Not(c))
        t0 = time.time()
        r = s.check()
        self.n_z3 += 1
        self.E.z3_seconds += time.time() - t0
        model = None
        if r == z3.sat:
            try:
                m = s.model()
                model = {str(d): m[d].as_long() for d in m.decls() if str(d).startswith(('v_', 's_'))}
            except Exception:
                model = None
        if r == z3.unsat:
            return True
        if self.derived and not must:
            self.deferred.append((cond, what))
            self.global_facts.append(c)
            return True
        if r == z3.sat:
            raise ProofFail(f'{what}: does not follow from the facts (counterexample found)', status='NOT_PROVED',
                            counterexample=model)
        raise ProofFail(f'{what}: Z3 returned {r} (time-out {self.E.z3_timeout_ms} ms)', status='UNKNOWN')

    def instantiate(self, e, env_var, env_ref=None, env_fold=None, memo=None):
        def fn(leaf):
            k = leaf[0]
            if k == 'var':
                return env_var[leaf[1]]
            if k == 'ref':
                return env_ref[leaf[1]]
            if k == 'fold':
                return env_fold(leaf)
            return None
        return X.subst(e, fn, memo)


class Engine:
    """Recursive prover over a Program (see the module docstring)."""

    def __init__(self, program, *, z3_timeout_ms=120_000, leaf_budget=600, leaf_max_mb=1500, log=None,
                 leaf_max_ops=400_000, seed=11, identity=None, leaf_prover=None):
        self.P = program
        self.z3_timeout_ms = z3_timeout_ms
        self.leaf_budget, self.leaf_max_mb, self.leaf_max_ops = leaf_budget, leaf_max_mb, leaf_max_ops
        self.log = log or (lambda m: None)
        self.seed = seed
        self.cache = {}
        self.records = []
        self.failed_leaves = []
        self.z3_seconds = 0.0
        self.identity = identity or _identity
        self.leaf_prover = leaf_prover or prove_leaf
        self.stats = dict(leaf_proofs=0, leaf_cached=0, compositions=0, compositions_cached=0, obligations=0, z3_calls=0)
        try:
            import z3  # noqa: F401
            self.z3 = z3
        except Exception:
            self.z3 = None

    # ---------------------------------------------------------------- entry
    def prove(self, node, contract=None, *, _depth=0):
        """-> dict(status, contract, mode, ...).  ``contract`` None: the adapter's registry, else derived."""
        if self.z3 is None:
            return dict(status='UNKNOWN', path=node.path, detail='z3 is not importable: no composition is proved')
        source = 'given'
        if contract is None:
            contract = self.P.adapter.contract(node)
            source = 'registry' if contract is not None else 'derived'
        try:
            qh, lh = self.identity(self.P, node)
        except Exception as e:
            return dict(status='ERROR', path=node.path, type=node.type,
                        detail=f'cannot cut the node: {type(e).__name__}: {e}')
        key = (qh, lh, contract.digest() if contract is not None else 'derived')
        hit = self.cache.get(key)
        if hit is not None:
            r = dict(hit)
            r.update(path=node.path, dedup_of=hit['path'])
            self.stats['leaf_cached' if hit.get('mode') == 'leaf' else 'compositions_cached'] += 1
            self.records.append({k: r.get(k) for k in ('path', 'type', 'status', 'mode', 'method', 'qasm_sha256', 'ops',
                                                        'contract_source', 'contract_name', 'detail', 'dedup_of')})
            return r
        children = None
        if contract is not None and contract.kind in ('primitive', 'box') and node.nops <= self.leaf_max_ops:
            mode = 'leaf'
        else:
            try:
                children = self.P.children(node)
            except Exception as e:  # noqa: BLE001
                return dict(status='UNKNOWN', path=node.path, type=node.type, mode='compose',
                            detail=f'children not available: {type(e).__name__}: {e}')
            mode = 'compose' if children else 'leaf'
        if mode == 'leaf' and contract is None:
            r = dict(status='UNKNOWN', mode='leaf', detail='no contract in the registry and no decomposition')
        elif mode == 'leaf':
            r = self._leaf(node, contract)
            if r['status'] not in ('PROVED', 'REFUTED'):
                try:
                    children = self.P.children(node)
                except Exception as e:  # noqa: BLE001  (e.g. the events below a primitive were not kept)
                    children = None
                    r['detail'] = (r.get('detail') or '') + f'; no fallback to the children: {e}'
                if children:
                    self.log(f'{node.path}: leaf proof {r["status"]} ({(r.get("detail") or "")[:80]}); '
                             f'trying the children')
                    r2 = self._compose(node, contract, children, _depth)
                    r2['leaf_attempt'] = {k: r.get(k) for k in ('status', 'detail', 'method')}
                    r = r2
        else:
            r = self._compose(node, contract, children, _depth)
        r.update(path=node.path, type=node.type, qasm_sha256=qh, layout_sha256=lh, contract_source=source,
                 ops=node.nops)
        if r['status'] == 'PROVED':
            self.cache[key] = r
        self.records.append({k: r.get(k) for k in ('path', 'type', 'status', 'mode', 'method', 'qasm_sha256', 'ops',
                                                    'contract_source', 'contract_name', 'detail')})
        return r

    # ---------------------------------------------------------------- leaves
    def _leaf(self, node, contract):
        self.stats['leaf_proofs'] += 1
        r = self.leaf_prover(self.P.node_ops(node), node.ins, node.outs, contract, budget=self.leaf_budget,
                             max_mb=self.leaf_max_mb, seed=self.seed, label=node.path)
        out = dict(status=r['status'], mode='leaf', method=r.get('method'), detail=r.get('detail'),
                   contract=contract, contract_name=contract.name, contract_sha256=contract.digest(),
                   sim_cases=r.get('sim_cases'), witness=r.get('witness'), leaf=r)
        if r['status'] != 'PROVED':
            self.failed_leaves.append((node, contract, r))
        self.log(f'leaf {node.path} [{contract.name}] {r["status"]} ops={r.get("ops")} {r.get("wall")}s')
        return out

    # ---------------------------------------------------------------- composition
    def _compose(self, node, contract, children, depth):
        self.stats['compositions'] += 1
        t0 = time.time()
        cp = _Comp(self, node, contract, self.z3)
        info = dict(proved_children=0, folds={})
        child_summary = {}
        derived = None
        try:
            if contract is not None:
                contract.validate()
                if node.widths_in() != contract.inputs or node.widths_out() != contract.outputs:
                    raise ProofFail(f'registers of the node {node.widths_in()} -> {node.widths_out()} do not match '
                                    f'the contract {contract.inputs} -> {contract.outputs}')
            self._coverage(node, children)
            for name, ws in node.ins.items():
                cp.write_value(ws, cp.new_val(('var', name), len(ws)))
            env_in = {name: ('var', name) for name in node.ins}
            node_wires = set(node.in_wires()) | set(node.out_wires()) | self.P.node_support(node)
            cp.touched |= set(node.in_wires())
            if contract is not None and contract.pre != 1:
                if contract.uses_ref_in_pre():
                    raise ProofFail('the precondition of a composite contract must read only its inputs', status='UNKNOWN')
                cp.global_facts.append(cp.T.cond(contract.pre))
            folds = {f.name: dict(f=f, i=0, emits={r: [] for r in f.emit})
                     for f in (contract.folds if contract is not None else [])}
            active = None
            for kid in children:
                if kid.kind == 'wires':
                    self._check_noop(kid)
                    if sorted(kid.in_wires()) != sorted(kid.out_wires()):
                        raise ProofFail(f'{kid.path}: a renaming node changes the set of wires')
                    continue
                if kid.kind == 'alloc':
                    self._check_noop(kid)
                    for q in kid.out_wires():
                        if q in cp.touched or q in cp.retired or cp.state.get(q, ZERO_C) != ZERO_C:
                            raise ProofFail(f'{kid.path}: allocated wire {q} is not fresh')
                        cp.touched.add(q)
                    continue
                if kid.kind == 'free':
                    self._check_noop(kid)
                    for name, ws in kid.ins.items():
                        for q in ws:
                            if q in cp.retired:
                                raise ProofFail(f'{kid.path}: wire {q} is freed twice')
                        cp.obligation(X.eq(cp.read(ws), 0), f'{kid.path}: freed register {name} is 0', must=True)
                        for q in ws:
                            cp.state[q] = ZERO_C
                            cp.retired.add(q)
                    continue
                # ---- a proper child (maybe an iteration of a loop)
                if active is not None and kid.type != folds[active]['f'].body:
                    if folds[active]['i'] < folds[active]['f'].count:
                        raise ProofFail(f'loop {active} interrupted by {kid.path} after {folds[active]["i"]} '
                                        f'of {folds[active]["f"].count} iterations')
                    active = None
                if active is None:
                    for fn, fs in folds.items():
                        if fs['i'] < fs['f'].count and kid.type == fs['f'].body:
                            active = fn
                            break
                fstate = folds[active] if active is not None else None
                if fstate is not None and fstate['i'] >= fstate['f'].count:
                    raise ProofFail(f'loop {active}: more iterations than its count {fstate["f"].count}')
                kc = fstate['f'].body_contract if fstate is not None else None
                kr = self.prove(kid, kc, _depth=depth + 1)
                sk = (kr.get('type'), kr.get('status'), kr.get('contract_name') or kr.get('contract_source'))
                child_summary[sk] = child_summary.get(sk, 0) + 1
                if kr['status'] != 'PROVED':
                    raise ProofFail(f'child {kid.path} ({kid.type}) is not proved: {kr.get("status")} '
                                    f'{(kr.get("detail") or "")[:300]}', status='NOT_PROVED', child=kr.get('path'),
                                    child_status=kr.get('status'))
                info['proved_children'] += 1
                if fstate is not None:
                    self._apply_iteration(cp, kid, fstate, env_in)
                else:
                    self._apply(cp, kid, kr['contract'])
            for fn, fs in folds.items():
                if fs['i'] != fs['f'].count:
                    raise ProofFail(f'loop {fn}: {fs["i"]} iterations found, the contract states {fs["f"].count}')
                info['folds'][fn] = dict(iterations=fs['i'], body=fs['f'].body, step_lemma=fs.get('lemma'))
            # ---- the end of the node: every non-output wire is 0
            out_wires = set(node.out_wires())
            for q in sorted(node_wires - out_wires):
                s = cp.state.get(q, ZERO_C)
                if s == ZERO_C:
                    continue
                if s[0] == 'c':
                    raise ProofFail(f'wire {q} ends in state 1 (it is not an output of the node)')
                cp.obligation(X.eq(X.bit(('val', s[0]), s[1]), 0), f'non-output wire {q} ends in 0', must=True)
            if contract is not None:
                env_fold = self._fold_env(cp, folds, env_in)
                memo, env_ref = {}, {}
                for nm, e in contract.defs:
                    env_ref[nm] = cp.instantiate(e, env_in, env_ref, env_fold, memo)
                for name, ws in node.outs.items():
                    want = cp.instantiate(contract.out[name], env_in, env_ref, env_fold, memo)
                    cp.obligation(X.eq(cp.read(ws), want), f'output {name} equals the contract', must=True)
            else:
                derived = self._derive(cp, node)
            status = 'PROVED'
            detail = (f'{info["proved_children"]} children proved, {cp.n_obl} obligations ({cp.n_z3} by Z3), '
                      f'{len(cp.deferred)} deferred to the parent')
            fail = None
        except ProofFail as e:
            status, detail, fail, derived = e.status, str(e), e.detail, None
        except (ContractError, X.ExprError) as e:
            status, detail, fail, derived = 'UNKNOWN', f'{type(e).__name__}: {e}', None, None
        self.stats['obligations'] += cp.n_obl
        self.stats['z3_calls'] += cp.n_z3
        used = contract if contract is not None else derived
        res = dict(status=status, mode='compose', method='z3_lia_composition', detail=detail, failure=fail,
                   contract=used, contract_name=used.name if used is not None else None,
                   contract_sha256=used.digest() if used is not None else None,
                   obligations=cp.n_obl, z3_calls=cp.n_z3, deferred=len(cp.deferred),
                   seconds=round(time.time() - t0, 2),
                   children=[dict(type=k[0], status=k[1], contract=k[2], count=v) for k, v in child_summary.items()],
                   folds=info.get('folds'))
        self.log(f'compose {node.path} {status} {detail[:200]}')
        return res

    def _check_noop(self, kid):
        if kid.span is not None and kid.span[1] > kid.span[0]:
            raise ProofFail(f'{kid.path}: a bookkeeping node owns operations')

    def _coverage(self, node, children):
        if node.info.get('contiguous') is False:
            raise ProofFail(f'coverage: the operations of {node.path} are not contiguous')
        if node.span is None:
            for kid in children:
                if kid.span is not None:
                    raise ProofFail(f'coverage: {kid.path} owns operations but its parent owns none')
            return
        a, b = node.span
        pos = a
        for kid in children:
            if kid.info.get('contiguous') is False:
                raise ProofFail(f'coverage: the operations of {kid.path} are not contiguous')
            if kid.span is None:
                continue
            ka, kb = kid.span
            if ka != pos:
                raise ProofFail(f'coverage: operations {pos}..{ka} of {node.path} belong to no child '
                                f'(or children overlap) before {kid.path}')
            if kb < ka:
                raise ProofFail(f'coverage: {kid.path} has a negative range')
            pos = kb
        if pos != b:
            raise ProofFail(f'coverage: operations {pos}..{b} of {node.path} belong to no child')

    def _fresh(self, cp, kid):
        ins = set(kid.in_wires())
        need_zero = (self.P.node_support(kid) | set(kid.out_wires())) - ins
        for q in need_zero:
            if q in cp.retired:
                raise ProofFail(f'freshness: {kid.path} uses wire {q} after it was freed')
            if cp.state.get(q, ZERO_C) != ZERO_C:
                raise ProofFail(f'freshness: {kid.path} assumes wire {q} is 0, but it carries a value '
                                f'(it is not an input of the child)')
        for q in ins:
            if q in cp.retired:
                raise ProofFail(f'wire flow: {kid.path} reads wire {q} after it was freed')
        return need_zero | ins

    def _apply(self, cp, kid, kc):
        if kc is None:
            raise ProofFail(f'{kid.path}: proved without a contract')
        if kid.widths_in() != kc.inputs or kid.widths_out() != kc.outputs:
            raise ProofFail(f'{kid.path}: registers do not match its contract')
        touched = self._fresh(cp, kid)
        env = {name: cp.read(ws) for name, ws in kid.ins.items()}
        memo, env_ref = {}, {}
        for nm, e in kc.defs:
            env_ref[nm] = cp.instantiate(e, env, env_ref, None, memo)
        if kc.pre != 1:
            pre = cp.instantiate(kc.pre, env, env_ref, None, memo)
            for j, conj in enumerate(_conjuncts(pre)):
                cp.obligation(conj, f'{kid.path}: precondition of {kc.name} (part {j})', must=False)
        self._write_outputs(cp, kid, kc, env, env_ref, memo, touched)

    def _write_outputs(self, cp, kid, kc, env, env_ref, memo, touched, override=None):
        snap = {q: cp.state.get(q, ZERO_C) for q in kid.in_wires()}
        new = {}
        ident = kc.identity_outputs()
        for name, ws in kid.outs.items():
            if override is not None and name in override:
                new[name] = ('vid', override[name])
            elif name in ident and len(kid.ins[name]) == len(ws):
                new[name] = ('copy', [snap[q] for q in kid.ins[name]])
            else:
                e = cp.instantiate(kc.out[name], env, env_ref, None, memo)
                new[name] = ('vid', cp.new_val(e, len(ws)))
        for q in touched:
            cp.state[q] = ZERO_C
        for name, ws in kid.outs.items():
            kind, v = new[name]
            if kind == 'copy':
                for q, s in zip(ws, v):
                    cp.state[q] = s
            else:
                cp.write_value(ws, v)
        cp.touched |= touched | set(kid.out_wires())
        return new

    # ---------------------------------------------------------------- loops
    def _apply_iteration(self, cp, kid, fs, env_in):
        f, i = fs['f'], fs['i']
        bc = f.body_contract
        if kid.widths_in() != bc.inputs or kid.widths_out() != bc.outputs:
            raise ProofFail(f'{kid.path}: registers do not match the body of loop {f.name}')
        touched = self._fresh(cp, kid)
        expected = {}
        for r in f.state:
            expected[r] = cp.instantiate(f.init[r], env_in) if i == 0 else ('sym', f'{f.name}.{r}.{i}')
        for r, e in f.inputs[i].items():
            expected[r] = cp.instantiate(e, env_in)
        for name, ws in kid.ins.items():
            cp.obligation(X.eq(cp.read(ws), expected[name]),
                          f'{kid.path}: loop {f.name} iteration {i}: input {name} is the stated one', must=True)
        memo, env_ref = {}, {}
        for nm, e in bc.defs:
            env_ref[nm] = cp.instantiate(e, expected, env_ref, None, memo)
        if i == 0:
            cp.obligation(cp.instantiate(f.invariant, expected), f'loop {f.name}: invariant before iteration 0',
                          must=True)
            fs['lemma'] = self._step_lemma(cp, f, expected)
        if bc.pre != 1:
            cp.obligation(cp.instantiate(bc.pre, expected, env_ref, None, memo),
                          f'{kid.path}: precondition of the loop body (iteration {i})', must=True)
        override, nxt = {}, {}
        for r in f.state:
            name = f'{f.name}.{r}.{i + 1}'
            cp.add_sym(name, bc.outputs[r])
            override[r] = cp.new_val(('sym', name), bc.outputs[r])
            nxt[r] = ('sym', name)
        env_next = dict(expected)
        env_next.update(nxt)
        cp.add_fact(cp.instantiate(f.invariant, env_next))
        new = self._write_outputs(cp, kid, bc, expected, env_ref, memo, touched, override)
        for r in f.emit:
            kind, v = new[r]
            if kind == 'copy':
                raise ProofFail(f'loop {f.name}: emitted register {r} is an unchanged input')
            fs['emits'][r].append(('val', v))
        fs['i'] = i + 1

    def _step_lemma(self, cp, f, expected0):
        """for every body input in range (fixed registers: their stated value): body pre => invariant after"""
        z3 = self.z3
        bc = f.body_contract
        fixed = f.fixed_registers()
        s = z3.Solver()
        s.set('timeout', self.z3_timeout_ms)
        env = {}
        for r, w in bc.inputs.items():
            if r in fixed:
                env[r] = expected0[r]
            else:
                nm = f'lemma.{f.name}.{r}'
                env[r] = ('sym', nm)
                cp.syms[nm] = w
                t = cp.T.term(env[r])
                s.add(t >= 0, t < (1 << w))
        # facts about the values that the fixed registers carry (their ranges, the parent facts)
        for r in fixed:
            s.add(*cp.facts_for(env[r]))
        for name, w in cp.node_widths.items():
            t = cp.T.term(('var', name))
            s.add(t >= 0, t < (1 << w))
        if cp.C is not None and cp.C.pre != 1:
            s.add(cp.T.cond(cp.C.pre))
        memo, env_ref = {}, {}
        for nm, e in bc.defs:
            env_ref[nm] = cp.instantiate(e, env, env_ref, None, memo)
        if bc.pre != 1:
            s.add(cp.T.cond(cp.instantiate(bc.pre, env, env_ref, None, memo)))
        nxt = dict(env)
        for r in f.state:
            nxt[r] = cp.instantiate(bc.out[r], env, env_ref, None, memo)
        s.add(z3.Not(cp.T.cond(cp.instantiate(f.invariant, nxt))))
        t0 = time.time()
        r = s.check()
        self.z3_seconds += time.time() - t0
        cp.n_z3 += 1
        cp.n_obl += 1
        if r != z3.unsat:
            raise ProofFail(f'loop {f.name}: step lemma (body precondition => invariant after the body) is not '
                            f'proved: Z3 {r}', status='NOT_PROVED' if r == z3.sat else 'UNKNOWN')
        return 'proved: Z3 unsat of the negation for a symbolic iteration'

    def _fold_env(self, cp, folds, env_in):
        def env(leaf):
            f, key, i = leaf[1], leaf[2], leaf[3]
            fs = folds[f]
            if key.startswith('emit:'):
                return fs['emits'][key[5:]][i]
            if i == 0:
                return cp.instantiate(fs['f'].init[key], env_in)
            return ('sym', f'{f}.{key}.{i}')
        return env

    # ---------------------------------------------------------------- derived contracts
    def _derive(self, cp, node):
        """the contract of a node without one: outputs = the state, pre = the deferred conditions"""
        defs, names = [], {}

        def lift(e):
            def fn(leaf):
                if leaf[0] == 'val':
                    vid = leaf[1]
                    v = cp.vals[vid]
                    if type(v.expr) is tuple and v.expr[0] == 'var':
                        return v.expr
                    if vid not in names:
                        inner = lift(v.expr)
                        nm = f'd{vid}'
                        defs.append((nm, inner))
                        names[vid] = ('ref', nm)
                    return names[vid]
                if leaf[0] == 'sym':
                    raise ProofFail('a derived contract cannot export a loop symbol', status='UNKNOWN')
                return None
            return X.subst(e, fn)
        out = {name: lift(cp.read(ws)) for name, ws in node.outs.items()}
        pre = 1
        if cp.deferred:
            dmap = {}
            conds = []
            for c, _ in cp.deferred:
                lc = lift(c)
                conds.append(lc)
            dmap = dict(defs)
            pre = X.all_(conds)
            pre = X.subst(pre, lambda leaf: _inline(leaf, dmap) if leaf[0] == 'ref' else None)
        c = Contract(inputs=node.widths_in(), outputs=node.widths_out(), out=out, pre=pre, defs=defs,
                     name=f'derived:{node.type}', kind='derived',
                     note='outputs = composition of the proved children; pre = their deferred preconditions')
        return c.validate()


def _conjuncts(e):
    out, stack = [], [e]
    while stack:
        x = stack.pop()
        if type(x) is tuple and x[0] == 'and':
            stack.append(x[2])
            stack.append(x[1])
        elif x != 1:
            out.append(x)
    return out


def _inline(leaf, dmap):
    e = dmap[leaf[1]]
    return X.subst(e, lambda l: _inline(l, dmap) if l[0] == 'ref' else None)


def _identity(P, node):
    """(QASM hash, layout hash) of the node program, renumbered as in 0.7 (montgomery._core.compact)"""
    from ..qualtran.montgomery._core import compact, make_layout
    from ..qualtran.bridge.emit import emit_qasm
    from ..qualtran.bridge.strict_proof import digest_json
    from hashlib import sha256
    nq, cops, wi, wo, _ = compact(list(P.node_ops(node)), dict(node.ins), dict(node.outs))
    return sha256(emit_qasm(nq, cops).encode()).hexdigest(), digest_json(make_layout(nq, wi, wo, cops))
