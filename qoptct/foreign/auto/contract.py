# SPDX-License-Identifier: Apache-2.0
"""Contracts (what a block does) and folds (loops with an invariant), as data.

``Contract``
    ``inputs``  {register: width} -- the registers the block reads (LEFT and THRU), unsigned integers;
    ``outputs`` {register: width} -- the registers it returns (RIGHT and THRU);
    ``defs``    [(name, expression)] let-bindings, used through ('ref', name);
    ``pre``     a 0/1 expression over the inputs and the definitions: the domain (default 1 = every input);
    ``out``     {register: expression}: the EXACT value of every output register; it must lie in [0, 2**width) on the
                domain (the leaf prover checks it; a value outside is a failed proof, never a silent truncation);
    ``folds``   loops of a composite block (``Fold``), referenced through ('fold', name, key, i) leaves.
    Every wire of the block that is not an output must be 0 at the end; this is part of every proof, not of the data.

``Fold``
    A run of ``count`` identical children ``body`` of a composite block.  ``body_contract`` is the contract every
    iteration is proved against; ``state`` are the body registers carried from one iteration to the next (output of
    iteration i = input of iteration i + 1, on the same physical wires); ``init`` their values before iteration 0;
    ``inputs[i]`` the expected values of the other body input registers at iteration i (expressions over the parent's
    inputs); ``emit`` the body output registers whose values are exported as ('fold', name, 'emit:<reg>', i);
    ``invariant`` an expression over the state registers and the registers whose input is the same in every iteration:
    it is proved by induction (base case at iteration 0, one step lemma for a symbolic iteration) and is the only thing
    known about the carried state.  Without a fold a loop is simply unrolled.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import expr as X

__all__ = ['Contract', 'Fold', 'ContractError']


class ContractError(ValueError):
    pass


def _ints(d, what):
    out = {}
    for k, v in dict(d).items():
        if not isinstance(k, str) or not k or type(v) is not int or v < 1:
            raise ContractError(f'{what}: register {k!r} needs a positive integer width')
        out[k] = v
    return out


@dataclass
class Fold:
    name: str
    body: str
    count: int
    body_contract: 'Contract'
    state: list
    init: dict
    inputs: list
    emit: list = field(default_factory=list)
    invariant: object = 1

    def validate(self, parent_inputs):
        bc = self.body_contract.validate()
        if type(self.count) is not int or self.count < 1 or len(self.inputs) != self.count:
            raise ContractError(f'fold {self.name}: count and inputs do not match')
        for r in self.state:
            if r not in bc.inputs or r not in bc.outputs or bc.inputs[r] != bc.outputs[r]:
                raise ContractError(f'fold {self.name}: state register {r} must be an input and an output of equal width')
        if set(self.init) != set(self.state):
            raise ContractError(f'fold {self.name}: init must give every state register')
        others = set(bc.inputs) - set(self.state)
        for e in self.init.values():
            X.validate(e, ('var',))
            if not X.refs_of(e, 'var') <= set(parent_inputs):
                raise ContractError(f'fold {self.name}: init reads unknown parent registers')
        for i, d in enumerate(self.inputs):
            if set(d) != others:
                raise ContractError(f'fold {self.name}: iteration {i} must give exactly the non-state inputs '
                                    f'{sorted(others)}')
            for e in d.values():
                X.validate(e, ('var',))
                if not X.refs_of(e, 'var') <= set(parent_inputs):
                    raise ContractError(f'fold {self.name}: iteration {i} reads unknown parent registers')
        for r in self.emit:
            if r not in bc.outputs:
                raise ContractError(f'fold {self.name}: emit register {r} is not a body output')
        X.validate(self.invariant, ('var',))
        fixed = self.fixed_registers()
        if not X.refs_of(self.invariant, 'var') <= set(self.state) | fixed:
            raise ContractError(f'fold {self.name}: the invariant may only read state registers and registers whose '
                                f'input is the same in every iteration ({sorted(fixed)})')
        return self

    def fixed_registers(self):
        if not self.inputs:
            return set()
        first = self.inputs[0]
        return {r for r in first if all(X.to_json(d[r]) == X.to_json(first[r]) for d in self.inputs[1:])}

    def to_dict(self):
        return dict(name=self.name, body=self.body, count=self.count, body_contract=self.body_contract.to_dict(),
                    state=list(self.state), init={k: X.to_json(v) for k, v in self.init.items()},
                    inputs=[{k: X.to_json(v) for k, v in d.items()} for d in self.inputs],
                    emit=list(self.emit), invariant=X.to_json(self.invariant))


@dataclass
class Contract:
    inputs: dict
    outputs: dict
    out: dict
    pre: object = 1
    defs: list = field(default_factory=list)
    folds: list = field(default_factory=list)
    name: str = ''
    kind: str = 'primitive'          # primitive | spec | body | derived | box
    note: str = ''

    def validate(self):
        if getattr(self, '_valid', False):
            return self
        self.inputs = _ints(self.inputs, 'inputs')
        self.outputs = _ints(self.outputs, 'outputs')
        if set(self.out) != set(self.outputs):
            raise ContractError(f'{self.name}: every output register needs exactly one expression (no unspecified junk)')
        allowed = ('var', 'ref', 'fold') if self.folds else ('var', 'ref')
        names = set()
        for f in self.folds:
            f.validate(self.inputs)
        fold_names = {f.name: f for f in self.folds}
        for nm, e in self.defs:
            if not isinstance(nm, str) or not nm or nm in names:
                raise ContractError(f'{self.name}: bad or repeated definition {nm!r}')
            X.validate(e, allowed)
            self._check_leaves(e, names, fold_names)
            names.add(nm)
        X.validate(self.pre, ('var', 'ref'))
        self._check_leaves(self.pre, names, {})
        for e in self.out.values():
            X.validate(e, allowed)
            self._check_leaves(e, names, fold_names)
        object.__setattr__(self, '_valid', True)
        return self

    def _check_leaves(self, e, names, folds):
        if not X.refs_of(e, 'var') <= set(self.inputs):
            raise ContractError(f'{self.name}: an expression reads unknown input registers')
        if not X.refs_of(e, 'ref') <= names:
            raise ContractError(f'{self.name}: an expression reads an undefined (or later) definition')
        for f, key, i in X.refs_of(e, 'fold'):
            fd = folds.get(f)
            if fd is None:
                raise ContractError(f'{self.name}: unknown fold {f}')
            if key.startswith('emit:'):
                if key[5:] not in fd.emit or not 0 <= i < fd.count:
                    raise ContractError(f'{self.name}: bad fold emit {key}[{i}]')
            elif key not in fd.state or not 0 <= i <= fd.count:
                raise ContractError(f'{self.name}: bad fold state {key}[{i}]')

    def uses_ref_in_pre(self):
        return bool(X.refs_of(self.pre, 'ref'))

    # ------------------------------------------------------------------ exact evaluation (independent checks)
    def in_domain(self, values):
        try:
            self.evaluate(values)
            return True
        except (ContractError, X.ExprError):
            return False

    def evaluate(self, values, *, check_pre=True):
        """values {input: int} -> {output: int}.  Raises ContractError outside the domain or for an output that
        does not fit its register."""
        self.validate()
        if set(values) != set(self.inputs):
            raise ContractError(f'{self.name}: inputs {sorted(values)} != {sorted(self.inputs)}')
        for k, v in values.items():
            if type(v) is not int or not 0 <= v < (1 << self.inputs[k]):
                raise ContractError(f'{self.name}: input {k}={v!r} outside its register')
        env_vals = dict(values)
        if check_pre and not self.uses_ref_in_pre():
            if X.evaluate(self.pre, lambda leaf: env_vals[leaf[1]]) != 1:
                raise ContractError(f'{self.name}: input outside the domain (precondition false)')
        folds = {}
        for f in self.folds:
            folds[f.name] = self._run_fold(f, env_vals)
        refs = {}

        def env(leaf):
            if leaf[0] == 'var':
                return env_vals[leaf[1]]
            if leaf[0] == 'ref':
                return refs[leaf[1]]
            f, key, i = leaf[1:]
            return folds[f][key][i]
        memo = {}
        for nm, e in self.defs:
            refs[nm] = X.evaluate(e, env, memo)
        if check_pre and self.uses_ref_in_pre():
            if X.evaluate(self.pre, env, memo) != 1:
                raise ContractError(f'{self.name}: input outside the domain (precondition false)')
        res = {k: X.evaluate(e, env, memo) for k, e in self.out.items()}
        for k, v in res.items():
            if not 0 <= v < (1 << self.outputs[k]):
                raise ContractError(f'{self.name}: output {k}={v} does not fit {self.outputs[k]} bits')
        return res

    def _run_fold(self, f, parent_vals):
        bc = f.body_contract
        st = {r: [X.evaluate(f.init[r], lambda leaf: parent_vals[leaf[1]])] for r in f.state}
        em = {'emit:' + r: [] for r in f.emit}
        for i in range(f.count):
            body_in = {r: st[r][i] for r in f.state}
            for r, e in f.inputs[i].items():
                body_in[r] = X.evaluate(e, lambda leaf: parent_vals[leaf[1]])
            inv = X.evaluate(f.invariant, lambda leaf: body_in[leaf[1]])
            if inv != 1:
                raise ContractError(f'fold {f.name}: invariant false before iteration {i}')
            body_out = bc.evaluate(body_in)
            for r in f.state:
                st[r].append(body_out[r])
            for r in f.emit:
                em['emit:' + r].append(body_out[r])
        st.update(em)
        return st

    # ------------------------------------------------------------------ identity
    def to_dict(self):
        self.validate()
        return dict(format='qoptct.auto.contract/1', name=self.name, kind=self.kind, inputs=self.inputs,
                    outputs=self.outputs, pre=X.to_json(self.pre),
                    defs=[[nm, X.to_json(e)] for nm, e in self.defs],
                    out={k: X.to_json(e) for k, e in self.out.items()},
                    folds=[f.to_dict() for f in self.folds], note=self.note,
                    junk_policy='every output explicit; every other wire of the block is 0 at the end')

    def digest(self):
        d = getattr(self, '_digest', None)
        if d is None:
            d = X.digest(self.to_dict())
            object.__setattr__(self, '_digest', d)
        return d

    def identity_outputs(self):
        """output registers whose expression is exactly the input register of the same name"""
        return {k for k, e in self.out.items() if e == ('var', k) and self.inputs.get(k) == self.outputs[k]}
