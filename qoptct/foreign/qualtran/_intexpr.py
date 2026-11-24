"""Read an integer written as a small arithmetic expression, such as ``2**255-19``.

Only integer literals, + - * ** << and parentheses are accepted; names, calls and attribute access are refused.
Exponents and results are bounded, so a short text cannot make the reader run for a long time.
"""
import ast

_MAX_BITS = 1 << 16       # results up to 65536 bits
_MAX_EXP = 1 << 16        # exponents and shifts up to 65536


def int_expr(text):
    """Return the integer value of ``text`` (an int is returned as is); raise ValueError otherwise."""
    if isinstance(text, int) and not isinstance(text, bool):
        return int(text)
    if not isinstance(text, str):
        raise ValueError(f'not an integer expression: {text!r}')
    if len(text) > 200:
        raise ValueError('integer expression too long')
    try:
        tree = ast.parse(text.strip(), mode='eval')
    except SyntaxError as e:
        raise ValueError(f'cannot read the number {text!r}') from e
    return _ev(tree.body, text)


def _check(v, text):
    if v.bit_length() > _MAX_BITS:
        raise ValueError(f'number too large in {text!r}')
    return v


def _ev(node, text):
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return _check(node.value, text)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        v = _ev(node.operand, text)
        return -v if isinstance(node.op, ast.USub) else v
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Pow, ast.LShift)):
        a, b = _ev(node.left, text), _ev(node.right, text)
        if isinstance(node.op, ast.Add):
            return _check(a + b, text)
        if isinstance(node.op, ast.Sub):
            return _check(a - b, text)
        if isinstance(node.op, ast.Mult):
            if a.bit_length() + b.bit_length() > _MAX_BITS + 1:
                raise ValueError(f'number too large in {text!r}')
            return _check(a * b, text)
        if b < 0 or b > _MAX_EXP:
            raise ValueError(f'exponent or shift out of range in {text!r}')
        if isinstance(node.op, ast.Pow):
            if abs(a) > 1 and (abs(a).bit_length() - 1) * b > _MAX_BITS:
                raise ValueError(f'number too large in {text!r}')
            return _check(a ** b, text)
        return _check(a << b, text)
    raise ValueError(f'only integers and + - * ** << ( ) are allowed: {text!r}')
