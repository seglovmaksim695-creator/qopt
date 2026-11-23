"""Qualtran bridge tests: the ported checks return evidence dictionaries; pytest must ignore the value."""
import pytest


@pytest.hookimpl(tryfirst=True)
def pytest_pyfunc_call(pyfuncitem):
    funcargs = pyfuncitem.funcargs
    args = {name: funcargs[name] for name in pyfuncitem._fixtureinfo.argnames}
    pyfuncitem.obj(**args)
    return True
