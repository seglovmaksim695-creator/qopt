"""Cache location: ``$QOPT_CACHE`` if set, else ``~/.cache/qopt`` (the same default as the rest of qoptct)."""
import os
from pathlib import Path


def cache_dir(*parts):
    """Return (and create) the cache directory, optionally a sub-directory ``parts``."""
    base = os.environ.get('QOPT_CACHE')
    if not base:
        base = os.path.join(os.path.expanduser('~'), '.cache', 'qopt')
    path = Path(base).joinpath(*parts)
    path.mkdir(parents=True, exist_ok=True)
    return path
