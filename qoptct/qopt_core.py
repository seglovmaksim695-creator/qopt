#!/usr/bin/env python3
"""qoptct 0.5.0 (project "qopt") -- verified Clifford+T circuit optimizer.
Install: pip install qoptct  (neural mode --net: pip install "qoptct[net]"). Python: import qoptct.

  qoptct optimize INPUT.qc [-o OUTPUT] [--mode cx|balanced|t] [--clean] [--meas] [--net] [--post] [--report r.json]
      Input: Feynman .qc (H, X, Z, T, T*, S, S*, P, P*, tof with 1-3 qubits, multi-qubit Z/Zd).
      --mode cx (default): unitary result with minimum CX (portfolio of rewriting programs + phase folding;
      --net adds neural block resynthesis). The result is ALWAYS proven equal to the input
      (<=10 qubits: full unitary; otherwise PyZX, then a numerical check in complex128 up to 26 qubits, then the
      exact path-sum prover of qoptct.verify);
      if equality is not proven the input circuit is returned.
      --mode balanced / t: measurement-based circuits (OpenQASM 3) with fewer T; every such output is checked
      by the independent simulator verify_qasm (<=22 qubits) or proven exactly (qoptct.verify) before it is written.
      Optional modules (v0.4): candidate sources qoptct.contrib.arith / bandc / gf2m (recognized arithmetic,
      logical-AND compiler for large circuits, GF(2^m) multipliers) and the CX post-pass qoptct.contrib.cxpass;
      every module result is checked again by the core, and a failing, slow or missing module is skipped
      (--no-arith, --no-bandc, --no-gf2m, --no-cxpass, --no-contrib, --contrib-sec, --budget).
      --prove: finally prove the written result exactly (qoptct.verify; MQT QCEC for unitary circuits if installed:
      pip install "qoptct[qcec]").
      --time-limit S (default 3600): all searches, including the unitary portfolio, stop at S seconds and the best
      result checked so far is written (in the worst case the input itself); only the final check of that result
      may follow. Output: OpenQASM 3 if OUTPUT ends in .qasm (also for a unitary result), otherwise .qc for a
      unitary result. Exit code 0: result written; 3: usage or input error (one "error: ..." line on stderr).
  qoptct prove CIRCUIT REF.qc [--clean] [--budget S]
      Exact equivalence proof of an optimized circuit (OpenQASM 3 or .qc) against its input; prints JSON, exit code
      0 PROVED, 1 REFUTED, 2 UNKNOWN, 3 usage or input error (missing file, malformed circuit).
  qoptct version
"""
import re
import argparse, importlib.util, json, os, signal, subprocess, sys, time


class _TO(BaseException):  # не Exception: программы эволюции глотают except Exception и висели часами
    pass


def _alarm(sig, frm):
    raise _TO()


_END = None  # absolute end of the run (--time-limit), set by main(); None = no overall limit
_LIMITS = []  # ends of the enclosing limited() calls (they may nest)
_PROOF_RESERVE = 0.0  # with --time-limit: seconds before the end kept free of portfolio programs for the proof


def time_left():
    """Seconds left until the overall time limit of the run (inf without one)."""
    return float("inf") if _END is None else _END - time.time()


def _within(seconds):
    """Seconds for one sub-stage: min(its own limit, time left of the run); None if no time is left (the sub-stage
    is then skipped). seconds <= 0 means no own limit (0 is returned when the run has no limit)."""
    left = time_left()
    if left <= 0.01:
        return None
    own = seconds if seconds and seconds > 0 else float("inf")
    s = min(own, left)
    return 0 if s == float("inf") else s


def limited(seconds, fn, *a):
    """Run fn(*a) for at most `seconds` (<= 0: no own limit), never beyond an enclosing limited() call and never
    beyond the overall time limit of the run (--time-limit); raise _TO otherwise. Calls may nest."""
    end = time.time() + seconds if seconds and seconds > 0 else float("inf")
    end = min([end] + _LIMITS[-1:] + ([_END] if _END is not None else []))
    old = signal.signal(signal.SIGALRM, _alarm)
    _LIMITS.append(end)
    if end != float("inf"):
        signal.setitimer(signal.ITIMER_REAL, max(0.001, end - time.time()))
    try:
        return fn(*a)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        _LIMITS.pop()
        signal.signal(signal.SIGALRM, old)
        if _LIMITS and _LIMITS[-1] != float("inf"):  # re-arm the enclosing limit (fires at once if it has passed)
            signal.setitimer(signal.ITIMER_REAL, max(0.001, _LIMITS[-1] - time.time()))
from pathlib import Path

__version__ = "0.5.0"
TORCH_MSG = ('the neural mode (--net / net=True) needs PyTorch, which is an optional dependency of qoptct.\n'
             'Install it with:  pip install "qoptct[net]"   (or: pip install torch)')


def require_torch():
    """Import torch for the neural mode or raise ImportError with installation instructions."""
    try:
        import torch  # noqa: F401
    except ImportError as e:
        raise ImportError(TORCH_MSG) from e


ROOT = Path(__file__).resolve().parent
EV = ROOT / "portfolio"
ELITE = ("c0001_combo2", "c0002_vphase_push", "c0005_phasepoly_graysynth", "c0006_combo2_permpair", "c0007_conj_rtof",
         "c0008_cla_collapse", "c0009_phasepoly_rtof4", "c0010_combo2_majcarry", "c0011_rtof3_vphase",
         "c0013_unify_permpair_vphase", "c0014_majcarry_ppresynth", "c0015_mcx_resynth", "c0017_rtof3_esop")
PROG_SEC = int(os.environ.get("QOPT_PROG_SEC", "120")); PROVE_SEC = int(os.environ.get("QOPT_PROVE_SEC", "600"))
FOLD_SEC = int(os.environ.get("QOPT_FOLD_SEC", "1800"))
PYZX_SEC = int(os.environ.get("QOPT_PYZX_SEC", "120"))  # базовой свёртке на больших схемах (hwb) нужно больше 120 с
ONE = {"H": "h", "X": "x", "T": "t", "T*": "tdg", "S": "s", "S*": "sdg", "P": "s", "P*": "sdg"}
OUT = {"h": "H", "x": "X", "t": "T", "tdg": "T*", "s": "S", "sdg": "S*"}


def load_qc(path):
    """Read a Feynman .qc circuit -> (wire names, gates). Malformed input raises ValueError (OSError if the file
    cannot be read) with a one-line message naming the file and line."""
    names, idx, gates = None, {}, []
    with open(path, encoding="utf-8") as f:
        try:
            lines = f.readlines()
        except UnicodeDecodeError as e:
            raise ValueError(f"{path}: not a text file ({e.reason})") from None
    for no, line in enumerate(lines, 1):
        p = line.split()
        if not p or p[0].startswith("#"):
            continue
        if p[0] == ".v":
            names = p[1:]; idx = {q: i for i, q in enumerate(names)}; continue
        if p[0].startswith(".") or p[0] in ("BEGIN", "END"):
            continue
        if names is None:
            raise ValueError(f"{path}, line {no}: gate before the .v line (not a .qc circuit?): {line.strip()[:80]!r}")
        bad = [q for q in p[1:] if q not in idx]
        if bad:
            raise ValueError(f"{path}, line {no}: wire {bad[0]!r} is not declared in the .v line: {line.strip()[:80]}")
        op, qs = p[0], tuple(idx[q] for q in p[1:])
        if not 1 <= len(set(qs)) <= 3 and op in ("Z", "Zd", "tof"):
            raise ValueError(f"{path}, line {no}: unsupported: {line.strip()} ({op} takes 1-3 wires)")
        if op in ("Z", "Zd"):
            qs = tuple(dict.fromkeys(qs))  # Z a b a = CZ a b (x*x = x); в cycle_17_3 такие есть
            gates.append(({1: "z", 2: "cz", 3: "ccz"}[len(qs)], qs))
        elif op == "tof":
            if len(set(qs)) != len(qs):
                raise ValueError(f"{path}, line {no}: tof with a repeated qubit: {line.strip()}")
            gates.append(({1: "x", 2: "cx", 3: "ccx"}[len(qs)], qs))
        elif op in ONE and len(qs) == 1:
            gates.append((ONE[op], qs))
        else:
            raise ValueError(f"{path}, line {no}: unsupported gate: {line.strip()}")
    if not names:
        raise ValueError(f"{path}: " + ("empty file" if not "".join(lines).strip() else
                                        "no .v line (not a .qc circuit?)" if names is None else
                                        "the .v line declares no wires"))
    return names, gates


def save_qc(path, names, circ):
    with open(path, "w") as f:
        f.write(".v " + " ".join(names) + "\n\nBEGIN\n")
        for g, q in circ:
            if g in ("z", "cz", "ccz"):  # .qc: Z with 1-3 qubits (read back by load_qc)
                f.write("Z " + " ".join(names[i] for i in q))
            else:
                f.write(("tof " + " ".join(names[i] for i in q)) if g == "cx" else (OUT[g] + " " + names[q[0]]))
            f.write("\n")
        f.write("END\n")


QASM_SUFFIXES = (".qasm", ".qasm3")


def write_unitary(path, names, circ):
    """Write a unitary result: OpenQASM 3 if `path` ends in .qasm/.qasm3 (readable by Qiskit and `qoptct prove`),
    otherwise .qc. -> "qasm" or "qc"."""
    if str(path).lower().endswith(QASM_SUFFIXES):
        Path(path).write_text(low_to_qasm(len(names), circ))
        return "qasm"
    save_qc(path, names, circ)
    return "qc"


def _mod(path):
    s = importlib.util.spec_from_file_location(Path(path).stem, path)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


def prove(c, ref, n):
    """True if equality is proven; otherwise False/None."""
    sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "proto2"))
    import bench_real as br
    if n <= 10:
        from qiskit.quantum_info import Operator
        import block_synth as bs
        return bool(Operator(bs.to_qc(c, n)).equiv(Operator(bs.to_qc(ref, n))))
    try:  # PyZX на больших схемах может идти часами — свой предел, потом численная проверка
        if limited(PYZX_SEC, lambda: br.ours_to_zx(c, n).verify_equality(br.ours_to_zx(ref, n))):
            return True
    except (Exception, _TO):
        pass
    if n <= 26:
        from verify import numeric_equal
        return bool(numeric_equal(c, ref, n))
    return None


def bench_sha_ok(path, name):
    """True if the input file is byte-identical to the benchmark `name` the hand-made c1 constructions were built for."""
    import hashlib
    man = ROOT / "contrib" / "c1" / "best_clean" / "MANIFEST.json"
    if not man.exists():
        return False
    want = json.loads(man.read_text()).get(name)
    return bool(want) and hashlib.sha256(Path(path).read_bytes()).hexdigest() == want


def ours_file(path, f):
    """Cached unitary start circuit (contrib/eng1/ours/<name>.json) for the gadget pipeline.
    The shipped cache is keyed by file name only, so it is used only if the input's SHA-256 matches;
    files in other directories (temporary ones created by run_sel) are returned unchanged."""
    import hashlib
    f = Path(f)
    if not f.exists():
        return None
    if f.parent.resolve() != (ROOT / "contrib" / "eng1" / "ours").resolve():
        return f
    want = json.loads(f.read_text()).get("qc_sha256")
    return f if want and want == hashlib.sha256(Path(path).read_bytes()).hexdigest() else None


def optimize(path, net=False, cxpass=False, budget=None, log=None):
    names, hl = load_qc(path)
    return optimize_hl(names, hl, net, cxpass, budget, log)


def optimize_hl(names, hl, net=False, cxpass=False, budget=None, log=None):
    """Same as optimize() for a Toffoli-level gate list (used by the Qiskit plugin).
    cxpass=True: also run the optional CX post-pass (qoptct.contrib.cxpass) on the proven result; its output is
    used only if (CX, T) drops, T does not grow and the equality is proven again (prove())."""
    names, ref, best, best_name = _optimize_hl(names, hl, net)
    if cxpass and best_name:
        c = cxpass_unitary(best, ref, len(names), budget or Budget(), log if log is not None else [])
        if c is not None:
            best, best_name = c, best_name + "+cxpass"
    return names, ref, best, best_name


def _optimize_hl(names, hl, net=False):
    if net:
        require_torch()
    sys.path.insert(0, str(EV))
    import evo_lib, bench_hook
    n = len(names)
    ref = evo_lib.lower_std(hl)
    # every sub-stage gets min(its own limit, time left of the run (--time-limit)); a sub-stage cut by the limit is
    # skipped; in the worst case the input itself is returned (best_name None)
    cands, skipped = [], []

    def prog_sec(own):  # the programs stop _PROOF_RESERVE seconds before the end: time for proving the best one
        sec = _within(own)
        if sec is None or _END is None or not _PROOF_RESERVE:
            return sec
        sec = min(sec or float("inf"), _END - _PROOF_RESERVE - time.time())
        return sec if sec > 0.01 else None

    sec = prog_sec(FOLD_SEC)
    if sec is None:
        skipped.append("input+fold")
    else:
        try:
            cands.append(("input+fold", limited(sec, evo_lib.fold_cancel, ref, n)))
        except _TO:
            pass
    for nm in ELITE:
        sec = prog_sec(PROG_SEC)
        if sec is None:
            skipped.append(nm[6:]); continue
        try:
            m = _mod(EV / "programs" / f"{nm}.py")
            cands.append((nm[6:], limited(sec, lambda: evo_lib.fold_cancel(
                [(g, tuple(q)) for g, q in bench_hook.run_start(m, hl, n)], n))))
        except (Exception, _TO) as e:
            print(f"  {nm}: error {e!r}"[:120], file=sys.stderr)
    if skipped:
        print(f"  portfolio: {len(skipped)} program(s) skipped (time limit)", file=sys.stderr)
    cands.sort(key=lambda x: evo_lib.key(x[1]))
    best_name, best = None, ref
    for name, c in cands[:3]:  # проверяем лучшие по очереди до первого доказанного
        sec = _within(PROVE_SEC)
        if sec is None:
            break
        t0 = time.time()
        try:
            ok = limited(sec, prove, c, ref, n)
        except _TO:
            ok = None
        if ok is None:  # out of reach of PyZX / the numerical check (> 26 qubits or slow): exact path-sum prover
            ok = _prove_unitary_exact(c, names, ref, None if sec == 0 else sec - (time.time() - t0))
        if ok is True:
            best_name, best = name, c; break
    sec = _within(0)
    if net and best_name and sec is not None:
        os.environ.update(GCANCEL="1", RTOF="1", PORTFOLIO="1", MODEL3="models/rt_3_si2_3q.pt")
        import bench_real as br
        br._init(str(ROOT / "models" / "mx_4_4q.pt"))
        try:
            r = limited(sec, br.job, ("q", n, best, False, 4))
            if evo_lib.key(r[3]) < evo_lib.key(best) and limited(_within(PROVE_SEC) or 0.001, prove, r[3], ref, n) is True:
                best_name, best = best_name + "+net", r[3]
        except _TO:
            pass
    return names, ref, best, best_name


def _prove_unitary_exact(c, names, ref, seconds):
    """Exact proof of a unitary candidate with the path-sum prover of qoptct.verify (when installed) -> True/None.
    Used when prove() gives no answer (PyZX failed within its limit and the circuit is too wide for the numerical
    check), so that the cx result on large inputs does not depend on the machine load."""
    if not contrib_available("verify"):
        return None
    sec = PROVE_SEC if seconds is None else seconds
    if _END is not None:
        sec = min(sec, time_left())
    if sec < 5:
        return None
    import tempfile, shutil
    tmp = Path(tempfile.mkdtemp(prefix="qopt_ref_"))
    try:
        path = str(tmp / "ref.qc")
        save_qc(path, names, ref)
        r, _ = call_contrib("verify", "prove", (low_to_qasm(len(names), c), path, False, float(sec)), seconds=sec)
        return True if isinstance(r, dict) and r.get("status") == "PROVED" else None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def optimize_meas(path):
    """Режим --meas: логическое И Гидни (вычисление 4 T, отмена измерением 0 T), места ищутся автоматически.
    -> (names, сегменты [(команды, MZ|None)], T, CX, доказано, как)."""
    import gidney, evo_lib
    names, hl0 = load_qc(path); n = len(names); ins = gidney.inputs_of(path)
    anc = set() if ins is None else {i for i, x in enumerate(names) if x not in ins}
    hl = evo_lib.fold_h(hl0)
    items, r1, r2 = gidney.transform(hl, n, anc)
    if not (r1 or r2):
        return names, None, None, None, None, "no AND sites"
    import meas_portfolio
    low = meas_portfolio.build_low(names, items, n)  # куски — портфелем qopt (каждый доказан), сегменты — свёрткой
    ok = gidney.anf_check(items, n, anc); how = "GF(2)+pieces proven"
    if ok is True and n <= 24:  # дополнительно — вся схема симуляцией с измерениями
        try:
            ok = limited(PROVE_SEC, gidney.verify, hl0, low, n, anc); how = "GF(2)+simulation"
        except _TO:
            pass
    if ok is None:
        try:
            ok = limited(PROVE_SEC, gidney.verify, hl0, low, n, anc); how = "simulation"
        except _TO:
            ok = None
    T = sum(evo_lib.tcount(s) for s, _ in low); CX = sum(evo_lib.cx(s) for s, _ in low)
    return names, low, T, CX, ok, how


QASM = {"h": "h", "x": "x", "z": "z", "s": "s", "sdg": "sdg", "t": "t", "tdg": "tdg", "cx": "cx", "cz": "cz"}


# ---- optional contributed modules (v0.4) ---------------------------------------------------------------------
# Candidate sources qoptct.contrib.{arith,bandc,gf2m}: candidates(qc_path, clean, budget_sec) -> [{"qasm", "note"}]
# CX post-pass qoptct.contrib.cxpass:           improve(qasm, clean, budget_sec) -> qasm
# Exact prover qoptct.verify:                   prove(qasm_or_path, ref_qc_path, clean, budget_sec) -> {"status", ...}
# Every module call runs in a child process with a hard wall-clock limit, so a crash, an exception or a hang of
# a module can never stop the optimizer; a missing module is skipped. Everything a module returns is checked
# again here (independent simulator up to SIM_MAX_QUBITS qubits, otherwise the exact prover) before it is used.
SOURCES = ("arith", "bandc", "gf2m")
CONTRIB_SEC = float(os.environ.get("QOPT_CONTRIB_SEC", "240"))
SIM_MAX_QUBITS = 22
# > SIM_MAX_QUBITS qubits: accept a candidate whose exact proof ends UNKNOWN (time/width) if the sampled checks of
# qoptct.verify pass; QOPT_REQUIRE_PROOF=1 or --require-proof accepts PROVED only
REQUIRE_PROOF = os.environ.get("QOPT_REQUIRE_PROOF", "0") == "1"
_GATE_RE = re.compile(r"(h|x|z|s|sdg|t|tdg|cx|cz) q\[\d+\](, q\[\d+\])*;")
_IF_RE = re.compile(r"if \(([\s()c\[\]\d&^]+)\) \{(.*)\}")  # AND (&&) and XOR (^) of measured bits
_COND_RE = re.compile(r"[\s()c\[\]\d&^!]+")  # conditions over measured bits (also negation, nested ifs)
_RUNNER = r"""
import json, sys, importlib
mod, fn, args_path, out_path = sys.argv[1:5]
a = json.load(open(args_path))
try:
    r = getattr(importlib.import_module(mod), fn)(*a["args"], **a["kwargs"])
    res = {"ok": True, "result": r}
except BaseException as e:
    import traceback
    res = {"ok": False, "error": repr(e)[:300], "trace": traceback.format_exc()[-1500:]}
with open(out_path, "w") as f:
    json.dump(res, f)
"""


class Budget:
    """Shared wall-clock budget of the optional module stages (None or <= 0: unlimited)."""

    def __init__(self, total=None):
        self.end = time.time() + total if total and total > 0 else None

    def give(self, want):  # also never beyond the overall time limit of the run (--time-limit)
        left = time_left() if self.end is None else min(self.end - time.time(), time_left())
        return max(0.0, min(want, left))


MODULES = {}  # optional override: short name ("contrib.arith", "verify", ...) -> importable module name


def _modname(name):
    return MODULES.get(name) or f"{__package__ or 'qoptct'}.{name}"


def contrib_available(name):
    """True if the optional module qoptct.<name> (e.g. "contrib.arith", "verify") is installed."""
    try:
        s = importlib.util.find_spec(_modname(name))
        return s is not None and s.origin is not None  # a bare directory (namespace package) does not count
    except (ImportError, ValueError):
        return False


def _kill_tree(pid):
    """SIGKILL all descendants of pid (also those that started their own session), best effort."""
    try:
        ps = subprocess.run(["ps", "-A", "-o", "pid=,ppid="], capture_output=True, text=True, timeout=10).stdout
    except Exception:  # noqa: BLE001
        return
    kids = {}
    for line in ps.splitlines():
        try:
            c, par = (int(x) for x in line.split())
        except ValueError:
            continue
        kids.setdefault(par, []).append(c)
    todo, seen = [pid], []
    while todo:
        for c in kids.get(todo.pop(), []):
            seen.append(c); todo.append(c)
    for c in seen:
        try:
            os.kill(c, signal.SIGKILL)
        except OSError:
            pass


def call_contrib(name, fn, args=(), kwargs=None, seconds=60.0, grace=None):
    """Run qoptct.<name>.<fn>(*args, **kwargs) in a child process, killed after seconds + grace
    (default grace: 5 s + 25 % of seconds, at most 30 s; covers interpreter start and imports).
    Returns (result, None) or (None, reason); never raises. Arguments and result must be JSON."""
    import tempfile
    if not contrib_available(name):
        return None, "not installed"
    if _END is not None:  # never beyond the overall time limit of the run (the final check runs without one)
        left = time_left()
        if left < 1:
            return None, "skipped (time limit)"
        seconds = min(seconds, left)
    grace = min(30.0, 5.0 + 0.25 * seconds) if grace is None else grace
    tmp = Path(tempfile.mkdtemp(prefix="qopt_mod_"))
    try:
        (tmp / "a.json").write_text(json.dumps({"args": list(args), "kwargs": kwargs or {}}))
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join([str(ROOT.parent)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
        p = subprocess.Popen([sys.executable, "-c", _RUNNER, _modname(name), fn,
                              str(tmp / "a.json"), str(tmp / "r.json")], env=env, cwd=tmp,
                             stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, start_new_session=True)
        try:
            _, err = p.communicate(timeout=seconds + grace)
        except BaseException as e:  # time-out, or an enclosing limited() / Ctrl-C: never leave the child running
            if not isinstance(e, subprocess.TimeoutExpired):
                _kill_tree(p.pid)
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except OSError:
                    p.kill()
                raise
            _kill_tree(p.pid)
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except OSError:
                p.kill()
            p.wait()
            return None, f"timeout after {seconds + grace:.0f} s"
        if not (tmp / "r.json").exists():
            return None, f"exit code {p.returncode}: " + err.decode(errors="replace")[-300:]
        r = json.loads((tmp / "r.json").read_text())
        return (r["result"], None) if r.get("ok") else (None, r.get("error", "error"))
    except Exception as e:  # noqa: BLE001 -- a module must never break the optimizer
        return None, repr(e)[:300]
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)


def qasm_stats(text):
    """(qubits, T, CX, measurements) of an OpenQASM 3 text in the subset written by qoptct, or None if it uses
    anything else (any other gate, register or statement). CX counts cx and cz, also inside conditionals."""
    N, T, CX, M = None, 0, 0, 0
    for L in text.splitlines():
        L = L.strip()
        if not L or L.startswith("//") or L.startswith("OPENQASM") or L.startswith("include"):
            continue
        m = re.fullmatch(r"qubit\[(\d+)\] q;", L)
        if m:
            if N is not None:
                return None
            N = int(m.group(1)); continue
        if re.fullmatch(r"bit\[\d+\] c;", L):
            continue
        if re.fullmatch(r"c\[\d+\] = measure q\[\d+\];", L):
            M += 1; continue
        if re.fullmatch(r"reset q\[\d+\];", L):
            continue
        body = [L]
        if L.startswith("if"):  # also nested ifs (the Qiskit-readable form of AND conditions)
            from qoptct.verify._parse import split_if
            try:
                parts = split_if(L)
            except ValueError:
                return None
            if not parts or any(not _COND_RE.fullmatch(c) for c, _ in parts):
                return None
            body = [g + ";" for _, gs in parts for g in gs]
        for b in body:
            if not _GATE_RE.fullmatch(b):
                return None
            g = b.split()[0]
            T += g in ("t", "tdg"); CX += g in ("cx", "cz")
    if N is None or max([int(x) for x in re.findall(r"q\[(\d+)\]", text)] or [0]) >= N:
        return None
    return N, T, CX, M


# ---- Qiskit-readable OpenQASM 3 -------------------------------------------------------------------------------
# qiskit.qasm3.loads (Qiskit 2.x) reads `if (c[k])` and `if (!c[k])`, also nested, but not the operators ^ and &&
# that the CX post-pass and the v0.3 compilers use in conditions. The written result is rewritten into that form
# (exactly the same operation on every branch) and proven again before it is used.
_SIMPLE_COND = re.compile(r"\(*\s*!?\s*c\[\d+\]\s*\)*")


def _simple_cond(c):
    c = c.strip()
    return bool(_SIMPLE_COND.fullmatch(c)) and c.count("(") == c.count(")")


def qiskit_readable(text):
    """True if every `if` condition (at every nesting level) is a single bit or its negation."""
    from qoptct.verify._parse import _if_parts, _IF_START
    for L in text.splitlines():
        s = L.strip()
        while True:
            i = s.find("if")
            while i >= 0 and not _IF_START.match(s[i:]):
                i = s.find("if", i + 1)
            if i < 0:
                break
            try:
                cond, block, rest = _if_parts(s[i:])
            except ValueError:
                return False
            if not _simple_cond(cond):
                return False
            s = block + " " + rest
    return True


_U1 = None


def _is_involution(gates):
    """True if the gate sequence (h x z s sdg t tdg cx cz) squares to the identity exactly."""
    import numpy as np
    global _U1
    if _U1 is None:
        r = 2 ** -0.5; w = np.exp(1j * np.pi / 4)
        _U1 = {"h": np.array([[r, r], [r, -r]]), "x": np.array([[0, 1], [1, 0]]), "z": np.diag([1, -1]),
               "s": np.diag([1, 1j]), "sdg": np.diag([1, -1j]), "t": np.diag([1, w]), "tdg": np.diag([1, np.conj(w)])}
    gs = []
    for g in gates:
        p = g.replace(",", " ").split()
        if not p or p[0] not in ("h", "x", "z", "s", "sdg", "t", "tdg", "cx", "cz"):
            return False
        gs.append((p[0], [int(x) for x in re.findall(r"q\[(\d+)\]", g)]))
    qs = sorted({q for _, qq in gs for q in qq})
    if len(qs) > 8:
        return False
    pos = {q: i for i, q in enumerate(qs)}; n = len(qs); U = np.eye(2 ** n, dtype=complex)
    for name, qq in gs:
        st = U.reshape([2] * n + [2 ** n])
        if name in _U1:
            st = np.moveaxis(np.tensordot(_U1[name], st, axes=([1], [pos[qq[0]]])), 0, pos[qq[0]])
        else:
            a, b = pos[qq[0]], pos[qq[1]]
            st = st.copy(); idx = [slice(None)] * (n + 1); idx[a] = 1
            sub = st[tuple(idx)]
            bb = b - (1 if b > a else 0)
            if name == "cx":
                sub = np.flip(sub, axis=bb)
            else:
                j = [slice(None)] * n; j[bb] = 1; sub = sub.copy(); sub[tuple(j)] *= -1
            st[tuple(idx)] = sub
        U = st.reshape(2 ** n, 2 ** n)
    return bool(np.allclose(U @ U, np.eye(2 ** n), atol=1e-9))


def _readable_if(cond, body):
    """Qiskit-readable statements equal to `if (cond) {body}` -> list of lines (minimum number of body copies)."""
    from qoptct.verify._parse import cond_poly, peval, bits
    if _simple_cond(cond):
        return ["if (%s) { %s }" % (cond.strip(), " ".join(g + ";" for g in body))], 1
    P = cond_poly(cond, lambda k: 1 << k)
    vs = sorted({k for m in P for k in bits(m)})
    opts = []
    if _is_involution(body):  # B^(m1 + m2 + ...) = B^(m1 xor m2 xor ...): one copy per ANF monomial
        opts.append([(["c[%d]" % k for k in bits(m)], body) for m in sorted(P)])
    if len(vs) <= 6:  # one copy per satisfying assignment (exactly one of them holds)
        sat = [a for a in range(1 << len(vs)) if peval(P, sum(1 << vs[i] for i in range(len(vs)) if a >> i & 1))]
        opts.append([(["c[%d]" % vs[i] if a >> i & 1 else "!c[%d]" % vs[i] for i in range(len(vs))], body) for a in sat])
    if not opts:
        raise ValueError("no Qiskit-readable form for the condition " + cond)
    best = min(opts, key=len)
    out = []
    for conds, b in best:
        if not conds:
            out.append(" ".join(g + ";" for g in b))
            continue
        t = "{ " + " ".join(g + ";" for g in b) + " }"
        for c in reversed(conds[1:]):
            t = "{ if (%s) %s }" % (c, t)
        out.append("if (%s) %s" % (conds[0], t))
    return out, len(best)


def _cost(gates):
    return (sum(g.split()[0] in ("t", "tdg") for g in gates), sum(g.split()[0] in ("cx", "cz") for g in gates))


def _parity_if(cond, body, t, bit, holds=None):
    """`if (cond) {body}` through a qubit q[t] in a known basis state: |0>, or |c[holds]> (measured into c[holds] and
    not touched since). X on q[t] for every ANF monomial of cond (xor c[holds]) turns it into |cond>; it is measured
    into the new bit c[bit]; the same X gates (or a reset correction) restore its state; then `if (c[bit]) {body}`.
    No T or CX; one more measurement."""
    from qoptct.verify._parse import cond_poly, bits
    P = cond_poly(cond, lambda k: 1 << k)
    if holds is not None:
        P = P ^ frozenset([1 << holds])

    def xs():
        r = []
        for m in sorted(P):
            ks = bits(m)
            if not ks:
                r.append("x q[%d];" % t); continue
            blk = "{ x q[%d]; }" % t
            for k in reversed(ks[1:]):
                blk = "{ if (c[%d]) %s }" % (k, blk)
            r.append("if (c[%d]) %s" % (ks[0], blk))
        return r
    restore = ["if (c[%d]) { x q[%d]; }" % (bit, t)] if holds is None else xs()
    return xs() + ["c[%d] = measure q[%d];" % (bit, t)] + restore + \
        ["if (c[%d]) { %s }" % (bit, " ".join(g + ";" for g in body))]


def make_qiskit_readable(text, parity=True):
    """Rewrite the conditions of a qoptct OpenQASM 3 text into the form read by qiskit.qasm3.loads (single bits,
    negation, nested ifs) -> new text with the same operation on every measurement branch, or None if impossible.
    Per condition the cheaper of two forms by (T, CX, measurements) is used: (a) copies of the body under nested
    single-bit ifs (one per ANF monomial if the body is an involution, else one per satisfying assignment);
    (b) the condition computed into a measured qubit that is known to be |0> (a measured qubit right after its
    reset correction, or one holding a just measured bit): no extra T/CX, one extra measurement (parity=False: only
    (a)). The result must be proven again before it is used."""
    from qoptct.verify._parse import split_if
    m = re.search(r"qubit\[(\d+)\]", text)
    mb = re.search(r"bit\[(\d+)\] c;", text)
    if not m:
        return None
    nbits = int(mb.group(1)) if mb else 0
    new_bits = 0
    zero, pending, out = set(), {}, []

    def touch(stmt):
        for q in re.findall(r"q\[(\d+)\]", stmt):
            zero.discard(int(q)); pending.pop(int(q), None)

    lines = text.splitlines()
    for i in range(len(lines) - 1):  # a reset correction right after an unreadable if moves before it (they commute:
        a_, b_ = lines[i].strip(), lines[i + 1].strip()  # other qubits), so its qubit is |0> for the parity form
        rb = re.fullmatch(r"if \(\(?c\[\d+\]\)?\) \{\s*x q\[(\d+)\];\s*\}", b_)
        if parity and rb and a_.startswith("if") and not qiskit_readable(a_) and \
                ("q[%s]" % rb.group(1)) not in a_:
            lines[i], lines[i + 1] = lines[i + 1], lines[i]
    try:
        for L in lines:
            s = L.strip()
            parts = split_if(s) if s.startswith("if") else None
            if parts is None:
                mm = re.fullmatch(r"c\[(\d+)\] = measure q\[(\d+)\];", s)
                touch(s)
                if mm:
                    pending[int(mm.group(2))] = int(mm.group(1))
                elif re.fullmatch(r"reset q\[(\d+)\];", s):
                    zero.add(int(re.search(r"\d+", s).group()))
                out.append(L); continue
            if len(parts) == 1 and _simple_cond(parts[0][0]) and qiskit_readable(s):
                cond, body = parts[0]
                rm = re.fullmatch(r"\(*\s*c\[(\d+)\]\s*\)*", cond.strip())
                tq = re.fullmatch(r"x q\[(\d+)\]", body[0].strip()) if len(body) == 1 else None
                is_reset = bool(rm and tq and pending.get(int(tq.group(1))) == int(rm.group(1)))
                touch(s)
                if is_reset:  # c[k] = measure q[t]; if (c[k]) { x q[t]; }  ->  q[t] is |0>
                    zero.add(int(tq.group(1)))
                out.append(L); continue
            for cond, body in parts:
                a_lines, copies = _readable_if(cond, body)
                ca = tuple(copies * x for x in _cost(body)) + (0,)
                cands = [(ca, a_lines, False)]
                bq = {int(q) for g in body for q in re.findall(r"q\[(\d+)\]", g)}
                free = sorted(zero - bq)
                held = sorted((q, k) for q, k in pending.items() if q not in bq)
                if parity and not _simple_cond(cond) and (free or held):
                    cb = _cost(body) + (1,)
                    par = _parity_if(cond, body, free[0], nbits + new_bits) if free else \
                        _parity_if(cond, body, held[0][0], nbits + new_bits, holds=held[0][1])
                    cands.append((cb, par, True))
                cost, lines, uses_bit = min(cands, key=lambda c: c[0])
                new_bits += uses_bit
                for q in {int(q) for g in body for q in re.findall(r"q\[(\d+)\]", g)}:
                    zero.discard(q); pending.pop(q, None)
                out += lines
    except ValueError:
        return None
    r = "\n".join(out) + "\n"
    if new_bits:
        if mb:
            r = r.replace(mb.group(0), "bit[%d] c;" % (nbits + new_bits), 1)
        else:
            r = r.replace(m.group(0) + " q;", m.group(0) + " q;\nbit[%d] c;" % new_bits, 1)
    return r if qiskit_readable(r) else None


def qasm_unitary(text):
    """(qubits, low-level gate list) of a measurement-free qasm text, else None."""
    st = qasm_stats(text)
    if st is None or st[3] or "if (" in text or "reset" in text:
        return None
    gates = []
    for L in text.splitlines():
        L = L.strip()
        if _GATE_RE.fullmatch(L):
            g, rest = L.split(None, 1)
            gates.append((g, tuple(int(x) for x in re.findall(r"q\[(\d+)\]", rest))))
    return st[0], gates


def low_to_qasm(n, c):
    return "\n".join(["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{n}] q;"] +
                     [f"{QASM[g]} " + ", ".join(f"q[{x}]" for x in q) + ";" for g, q in c]) + "\n"


def check_qasm(text, path, clean, budget=None, log=None):
    """Independent check of a candidate against the input .qc: simulator verify_qasm for <= SIM_MAX_QUBITS
    qubits (strict model, or the clean-ancilla model with clean=True), otherwise the exact prover
    (qoptct.verify, PROVED only). -> (ok, how)."""
    import tempfile, verify_qasm
    st = qasm_stats(text)
    if st is None:
        return False, "unsupported qasm"
    n = len(load_qc(path)[0])
    if st[0] < n:
        return False, "fewer qubits than the input"
    if st[0] <= SIM_MAX_QUBITS and ("^" in text or "reset" in text):
        # XOR conditions / resets (CX post-pass output): read by the dense checker of qoptct.verify, not by verify_qasm
        q, err = call_contrib("verify", "check_dense", (text, str(Path(path).resolve())), {"clean": bool(clean)},
                              seconds=PROVE_SEC)
        if err or not isinstance(q, dict) or q.get("F") is None:
            return False, f"dense check: {err or (q or {}).get('error')}"[:200]
        if not (q.get("ok") and q["F"] > 1 - 1e-9):
            return False, f"dense check (qoptct.verify) F={q['F']:.12f}"
        return _also_prove(text, path, clean, budget, f"dense check F={q['F']:.12f}")
    if st[0] <= SIM_MAX_QUBITS:
        with tempfile.NamedTemporaryFile("w", suffix=".qasm", delete=False) as f:
            f.write(text)
        try:
            F = limited(PROVE_SEC, verify_qasm.check, f.name, path, 30, 1, True, clean)
        except (Exception, _TO) as e:
            return False, f"simulator: {e!r}"[:200]
        finally:
            os.unlink(f.name)
        if not F > 1 - 1e-9:
            return False, f"simulator F={F:.12f}"
        # the simulator samples measurement outcomes; an exact proof (qoptct.verify) also covers every branch
        return _also_prove(text, path, clean, budget, f"simulator F={F:.12f}")
    sec = (budget or Budget()).give(PROVE_SEC)
    if sec < 5:
        return False, "no time left for the exact proof"
    ref = str(Path(path).resolve())
    r, err = call_contrib("verify", "prove", (text, ref, clean, sec), seconds=sec)
    if err or not isinstance(r, dict):
        return False, f"prover: {err or 'bad result'}"
    if r.get("status") == "PROVED":
        _PROVEN[_text_key(text, ref, clean)] = r.get("method")
        return True, f"prover PROVED ({r.get('method')})"
    if r.get("status") == "UNKNOWN" and not REQUIRE_PROOF:
        # proof out of reach (time/width): accept on the independent sampled checks of qoptct.verify
        # (sparse simulator / path-sum sampler, measurement outcomes sampled); the report says it is not a proof
        ok, how = strong_check(text, ref, clean)
        return ok, "prover UNKNOWN; " + how
    return False, f"prover {r.get('status')} ({r.get('method')})"


MODULE_SEC = {}  # per-module time limits in seconds (--source-sec NAME=SEC); default: module_sec()


def module_sec(name, path=None):
    """Time limit of one call of module `name`: --source-sec, else CONTRIB_SEC; the logical-AND compiler bandc gets
    more on large inputs (compilation + sampled checks of 1000+ wires take 6-10 minutes)."""
    if name in MODULE_SEC:
        return MODULE_SEC[name]
    if name == "bandc" and path is not None:
        try:
            names, hl = load_qc(path)
            if len(names) > 200 or len(hl) > 2000:
                return max(CONTRIB_SEC, 600.0)
            if len(names) > 22 or len(hl) > 500:
                return max(CONTRIB_SEC, 300.0)
        except Exception:  # noqa: BLE001
            pass
    return CONTRIB_SEC


SMALL_PROVE_SEC = float(os.environ.get("QOPT_SMALL_PROVE_SEC", "60"))
STRONG_SEC = float(os.environ.get("QOPT_STRONG_SEC", "900"))
_STRONG = {}  # sha256 of (text, ref, clean) -> (ok, how): a text is checked strongly only once per run
_PROVEN = {}  # sha256 of (text, ref, clean) -> prover method: texts proven exactly (PROVED) earlier in this process


def _text_key(text, ref, clean):
    import hashlib
    return hashlib.sha256(f"{text}\0{Path(ref).resolve()}\0{bool(clean)}".encode()).hexdigest()


def strong_trials(n_qubits):
    """Random input states (each with freshly sampled measurement outcomes) of the strong sampled check."""
    return 256 if n_qubits <= 14 else 128 if n_qubits <= 18 else 64 if n_qubits <= SIM_MAX_QUBITS else 32


def strong_check(text, ref, clean, seconds=None):
    """Sampled check with many measurement branches, used when the exact proof of qoptct.verify ends UNKNOWN
    (not a proof; the report says so). <= SIM_MAX_QUBITS qubits: dense state vector (qoptct.verify.check_dense,
    reads XOR conditions and resets) over strong_trials() random inputs; larger: sparse branch simulator, then (strict
    model) the sampled path-sum checker. A dropped conditional correction on a measured bit survives one trial with
    probability 1/2, so it survives the check with probability <= 2^-32. -> (ok, how)."""
    import hashlib
    key = hashlib.sha256(f"{text}\0{Path(ref).resolve()}\0{bool(clean)}".encode()).hexdigest()
    if key in _STRONG:
        return _STRONG[key]
    m = re.search(r"qubit\[(\d+)\]", text)  # the checkers of qoptct.verify read a wider subset than qasm_stats
    if m is None:
        return False, "strong check: no qubit register"
    nq = int(m.group(1))
    sec = STRONG_SEC if seconds is None else seconds
    ref = str(Path(ref).resolve())
    nt = strong_trials(nq)
    if nq <= SIM_MAX_QUBITS:
        q, err = call_contrib("verify", "check_dense", (text, ref), {"clean": bool(clean), "trials": nt, "seed": 7919},
                              seconds=sec)
        eng = f"dense {nt} trials"
    else:
        q, err = call_contrib("verify", "check_sparse", (text, ref), {"clean": bool(clean), "trials": nt, "seed": 7919},
                              seconds=sec)
        eng = f"sparse {nt} trials"
        if not err and isinstance(q, dict) and q.get("F") is None and not clean:
            q, err = call_contrib("verify", "check_big", (text, ref), {"x_trials": 8, "m_trials": 8, "seed": 7919},
                                  seconds=sec)
            eng = "path-sum sampler 8x8"
    if not err and isinstance(q, dict) and q.get("ok") and q.get("F") is not None and q["F"] > 1 - 1e-9:
        res = True, f"sampled check {eng} F={q['F']:.12f} (not a proof)"
    else:
        res = False, f"sampled check {eng} failed ({err or (q or {}).get('F', (q or {}).get('error'))})"[:300]
    if not err:  # a check stopped by a time limit is not a verdict: the final check of the run may repeat it
        _STRONG[key] = res
    return res


def _also_prove(text, path, clean, budget, how):
    """After a passed simulator check (<= SIM_MAX_QUBITS qubits): exact proof with a short budget when qoptct.verify
    is installed. REFUTED rejects the candidate; UNKNOWN (or no time for the proof) needs the strong sampled check
    (many measurement branches), since the 30-sample simulator alone can miss a dropped conditional correction."""
    if not contrib_available("verify"):
        return True, how
    sec = (budget or Budget()).give(SMALL_PROVE_SEC)
    st = None
    if sec >= 5:
        r, err = call_contrib("verify", "prove", (text, str(Path(path).resolve()), bool(clean), sec), seconds=sec)
        st = r.get("status") if isinstance(r, dict) else None
        if st == "REFUTED":
            return False, how + f"; prover REFUTED ({r.get('method')})"
        if st == "PROVED":
            _PROVEN[_text_key(text, path, clean)] = r.get("method")
            return True, how + f"; prover PROVED ({r.get('method')})"
    ok, h2 = strong_check(text, path, clean)
    return ok, how + f"; prover {st or 'not run'}; " + h2


def contrib_candidates(path, clean, sources, budget, log):
    """Verified candidates of the source modules: [((T, CX, extra_qubits), qasm, "module: note")]."""
    out, n = [], len(load_qc(path)[0])
    for name in sources:
        sec = budget.give(module_sec(name, path))
        if sec < min(5.0, CONTRIB_SEC):
            log.append(f"{name}: skipped (time budget)"); continue
        t0 = time.time()
        r, err = call_contrib(f"contrib.{name}", "candidates", (str(Path(path).resolve()), bool(clean)),
                              {"budget_sec": sec}, seconds=sec)
        if err:
            log.append(f"{name}: {err}"); continue
        good = 0
        cl = [(c, c.get("qasm") if isinstance(c, dict) else None) for c in (r if isinstance(r, list) else [])[:50]]
        cl = [(c, txt, qasm_stats(txt) if isinstance(txt, str) else None) for c, txt in cl]
        cl.sort(key=lambda x: (x[2] is None, x[2][1:3] if x[2] else ()))  # best (T, CX) first
        for c, txt, st in cl:
            if st is None:
                log.append(f"{name}: rejected candidate (unsupported qasm)"); continue
            if time_left() < 1:
                log.append(f"{name}: remaining candidates not checked (time limit)"); break
            ok, how = check_qasm(txt, path, clean, budget)
            if not ok:
                log.append(f"{name}: rejected candidate T={st[1]} CX={st[2]} ({how})"); continue
            good += 1
            out.append(((st[1], st[2], st[0] - n), txt, f"{name}: {str(c.get('note', ''))[:120]} [{how}]"))
        log.append(f"{name}: {good} verified of {len(r) if isinstance(r, list) else 0} in {time.time() - t0:.1f} s")
    return out


def cxpass_qasm(text, path, clean, budget, log, better):
    """Run qoptct.contrib.cxpass.improve on a measurement-based circuit; returns the improved text only if
    better(old_stats, new_stats) and it passes check_qasm, else None."""
    st0 = qasm_stats(text)
    sec = budget.give(module_sec("cxpass"))
    if st0 is None or sec < min(5.0, CONTRIB_SEC):
        return None
    r, err = call_contrib("contrib.cxpass", "improve", (text, bool(clean)),
                          {"budget_sec": sec, "ref": str(Path(path).resolve())}, seconds=sec)
    if err or not isinstance(r, str):
        log.append(f"cxpass: {err or 'bad result'}"); return None
    st = qasm_stats(r)
    if st is None or st[0] > st0[0] or not better(st0, st):
        log.append(f"cxpass: no improvement ({st0[1]}/{st0[2]} -> {st and st[1]}/{st and st[2]})"); return None
    ok, how = check_qasm(r, path, clean, budget)
    log.append(f"cxpass: {st0[1]}/{st0[2]} -> {st[1]}/{st[2]} T/CX, {'accepted' if ok else 'rejected'} ({how})")
    return r if ok else None


def cxpass_unitary(best, ref, n, budget, log, path=None):
    """CX post-pass on a unitary result; accepted only if (CX, T) drops, T does not grow and the equality is proven
    again (prove(), else the exact prover qoptct.verify). `path`: the input .qc (written from `ref` if not given),
    handed to the post-pass as its reference."""
    import evo_lib, tempfile, shutil
    if not contrib_available("contrib.cxpass"):
        return None
    sec = budget.give(module_sec("cxpass"))
    if sec < min(5.0, CONTRIB_SEC):
        return None
    tmp = None
    try:
        if path is None:
            tmp = Path(tempfile.mkdtemp(prefix="qopt_ref_"))
            path = str(tmp / "ref.qc")
            save_qc(path, [f"q{i}" for i in range(n)], ref)
        path = str(Path(path).resolve())
        r, err = call_contrib("contrib.cxpass", "improve", (low_to_qasm(n, best), False), {"budget_sec": sec, "ref": path},
                              seconds=sec)
        u = qasm_unitary(r) if isinstance(r, str) else None
        if err or u is None or u[0] != n:
            log.append(f"cxpass: {err or 'no unitary result'}"); return None
        c = evo_lib.lower_std(u[1])  # z -> s s, cz -> h cx h: the .qc writer and the prover know only h x s sdg t tdg cx
        if not (evo_lib.key(c) < evo_lib.key(best) and evo_lib.tcount(c) <= evo_lib.tcount(best)):
            log.append("cxpass: no improvement"); return None
        try:
            ok = limited(PROVE_SEC, prove, c, ref, n)
        except _TO:
            ok = None
        how = "proven"
        if ok is not True and contrib_available("verify"):
            ps = float(min(PROVE_SEC, max(1.0, time_left())))
            pr, _ = call_contrib("verify", "prove", (low_to_qasm(n, c), path, False, ps), seconds=ps)
            ok = isinstance(pr, dict) and pr.get("status") == "PROVED"
            how = f"qoptct.verify: {pr.get('status') if isinstance(pr, dict) else 'error'}"
        log.append(f"cxpass: CX {evo_lib.cx(best)} -> {evo_lib.cx(c)}, {how if ok is True else 'not proven (' + how + ')'}")
        return c if ok is True else None
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)


def save_qasm3(path, n, segs):
    """Circuit with measurements as OpenQASM 3 (measurement + conditional CZ correction + reset)."""
    m = sum(1 for _, mz in segs if mz)
    L = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{n}] q;"] + ([f"bit[{m}] c;"] if m else [])
    k = 0
    for seg, mz in segs:
        for g, q in seg:
            L.append(f"{QASM[g]} " + ", ".join(f"q[{x}]" for x in q) + ";")
        if mz:
            a, b, t = mz
            L += [f"c[{k}] = measure q[{t}];", f"if (c[{k}]) {{ cz q[{a}], q[{b}]; x q[{t}]; }}"]; k += 1
    Path(path).write_text("\n".join(L) + "\n")


def run_sel(path, names, best, cap, seconds, qout):
    """Выборочная гаджетизация H (агент tg2): минимум T при CX <= cap. Возвращает dict с T, CX, extra_qubits, fidelity
    (независимая проверка verify_qasm) или None."""
    import hashlib, tempfile, shutil, verify_qasm
    spec = importlib.util.spec_from_file_location("tg2_sel", ROOT / "contrib" / "tg2" / "sel.py")
    sel = importlib.util.module_from_spec(spec); spec.loader.exec_module(sel)
    h = hashlib.sha1((Path(path).read_text() + repr(best)).encode()).hexdigest()[:10]
    tmp = Path(tempfile.mkdtemp(prefix="qsel_")); stem = Path(path).stem + "_" + h
    qc = tmp / (stem + ".qc"); shutil.copy(path, qc); (tmp / (stem + ".json")).write_text("{}")
    ours = [(g, tuple(q)) for g, q in best]
    def get_start(name, start):
        nm, hl = load_qc(str(qc))
        return (len(nm), ours, None) if start == "ours" else (len(nm), sel.prep_gates(hl), None)
    sel.get_start = get_start; sel.OURS = tmp; sel.verify = lambda *x, **k: None
    try:
        r = sel.improve(str(qc), qasm_out=qout, cx_cap=cap, log=lambda *x: None, deadline=time.time() + seconds)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if r.get("extra_qubits", 99) + len(names) > 22:
        return None
    try:
        r["fidelity"] = limited(PROVE_SEC, verify_qasm.check, qout, path)
    except _TO:
        return None
    return r


AND_SEC = float(os.environ.get("QOPT_AND_SEC", "1800"))  # --and-sec: time of the v0.3 AND-compiler variant search
TIME_LIMIT = float(os.environ.get("QOPT_TIME_LIMIT", "3600"))  # --time-limit: overall seconds of an optimize run


def _stage(log, name, fn, *a, default=None):
    """Run one search stage of main() within the overall time limit; a stage cut by the limit returns `default`
    (its partial work is dropped). Other exceptions propagate as before."""
    left = time_left()
    if left < 1:
        log.append(f"{name}: skipped (time limit)"); return default
    try:
        return limited(0 if left == float("inf") else left, fn, *a)
    except _TO:
        log.append(f"{name}: stopped at the time limit"); return default


def _and_legacy(path, names, cap, variants, clean):
    """Candidates of the v0.3 logical-AND compilers: ({(T, CX, extra): qasm}, ok_keys, n)."""
    sys.path.insert(0, str(ROOT / "contrib" / "c2"))
    spec = importlib.util.spec_from_file_location("c2_run", ROOT / "contrib" / "c2" / "run.py")
    c2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c2)
    try:  # the variant search of c2 has no natural end on large inputs (ham15-low: hours): AND_SEC limits it
        cands, n = c2.candidates(path, variants, None, clean=clean,
                                 deadline=time.time() + AND_SEC if AND_SEC > 0 else None)
    except Exception:
        cands, n = None, len(names)
    # recognized multi-controlled Toffolis / ESOP (agent c12): AND ladder, RCCX rewrite, Qiskit KG24 candidates
    try:
        sys.path.insert(0, str(ROOT / "contrib" / "c12"))
        spec12 = importlib.util.spec_from_file_location("c12_mcx", ROOT / "contrib" / "c12" / "mcx.py")
        m12 = importlib.util.module_from_spec(spec12); spec12.loader.exec_module(m12)
        extra = m12.candidates(path, clean=clean)
        if extra:
            cands = dict(cands or {})
            for k12, q12 in extra.items():
                cands[tuple(k12[:3])] = q12
    except Exception:
        pass
    cands = dict(cands or {})
    ok = [k for k in cands if k[1] <= cap and n + k[2] <= SIM_MAX_QUBITS]
    if not ok:
        # circuits with non-classical parts: hybrid of agent c8 (classical blocks -> AND, H/phases stay)
        sys.path[:0] = [str(ROOT / "contrib" / "c8")]
        import hyb, andc
        n, parts, front = hyb.compile_hybrid(path, max(4, variants // 10), log=lambda *x: None)
        best = None
        for (_, _, _, ch) in front:
            txt, N = hyb.emit(n, parts, ch); T2, CX2 = andc.count(txt)
            if CX2 <= cap and N <= SIM_MAX_QUBITS and (best is None or (T2, CX2) < best[:2]):
                best = (T2, CX2, N - n, txt)
        if best is not None:
            cands = {best[:3]: best[3]}; ok = [best[:3]]
    if clean:  # hand-made constructions of agent c1 for specific benchmarks (only for byte-identical inputs)
        bc = ROOT / "contrib" / "c1" / "best_clean" / (Path(path).stem + ".qasm")
        if bc.exists() and bench_sha_ok(path, Path(path).stem):
            t = bc.read_text(); N = int(re.search(r"qubit\[(\d+)\]", t).group(1))
            k1 = (len(re.findall(r"\b(t|tdg)\s+q", t)), len(re.findall(r"\b(cx|cz)\s+q", t)), N - n)
            if k1[1] <= cap and N <= SIM_MAX_QUBITS:
                cands[k1] = t; ok = list(ok) + [k1]
    return cands, ok, n


def _and_finish(path, text, key, qout, clean, fidelity=None):
    """Write `text` to qout, check it (unless already checked: fidelity given) and polish CX with agent c1
    (strict model only; T does not grow, only checked results are taken). -> result dict or None."""
    import verify_qasm
    Path(qout).write_text(text)
    if fidelity is None:
        try:
            fidelity = float(limited(PROVE_SEC, lambda: verify_qasm.check(qout, path, clean=clean)))
        except _TO:
            return None
    res = dict(T=key[0], CX=key[1], extra_qubits=key[2], fidelity=float(fidelity))
    try:
        if clean:
            raise RuntimeError("c1 polish verifies in the strict model only")
        if key[2] + len(load_qc(path)[0]) > SIM_MAX_QUBITS:
            raise RuntimeError("c1 polish needs the simulator")
        sys.path.insert(0, str(ROOT / "contrib" / "c1"))
        spec1 = importlib.util.spec_from_file_location("c1_run", ROOT / "contrib" / "c1" / "run.py")
        c1 = importlib.util.module_from_spec(spec1); spec1.loader.exec_module(c1)
        q1 = qout + ".c1.qasm"
        r1 = c1.improve(path, q1, mode="c2", src=qout)
        if r1.get("fidelity") and r1["fidelity"] > 1 - 1e-9 and r1["T"] <= res["T"] and r1["CX"] < res["CX"]:
            import shutil; shutil.copy(q1, qout)
            res.update(T=r1["T"], CX=r1["CX"], fidelity=float(r1["fidelity"]), post="c1")
    except Exception:
        pass
    return res


def run_and(path, names, cap, variants, qout, clean=False, sources=None, budget=None, log=None, mods=None):
    """Logical-AND compilers: Toffolis via Gidney's logical AND into fresh ancillas + measurement-based
    uncomputation (H gates are not gadgetized). Minimum T with CX <= cap.
    v0.3 compilers (agents c2, c12, c8, c1; only circuits <= 22 qubits, independent simulator check) plus, from
    v0.4, the optional source modules `sources` (subset of SOURCES; every candidate is checked by check_qasm).
    The v0.3 result is computed exactly as before, so adding modules can only improve the result."""
    log = log if log is not None else []
    budget = budget or Budget()
    # the modules run first: if the overall time limit cuts the (open-ended) v0.3 search, their results are kept
    # (main() computes them before all v0.3 measurement stages and hands them in as `mods`)
    if mods is None:
        mods = contrib_candidates(path, clean, sources or (), budget, log) if sources else []
    legacy, ok, n = ({}, [], len(names))
    if variants > 0:
        try:
            legacy, ok, n = _stage(log, "and", _and_legacy, path, names, cap, variants, clean, default=({}, [], len(names)))
        except Exception as e:  # noqa: BLE001
            log.append(f"and: {e!r}"[:200])
    results = []
    for k in sorted(ok, key=lambda k: (k[0], k[1], k[2]))[:1]:  # v0.3: the minimum is checked, no fallback
        r = _stage(log, "and: check", _and_finish, path, legacy[k], k, qout + ".legacy.qasm", clean)
        if r:
            r["source"] = "and"; results.append((r, qout + ".legacy.qasm"))
    mods = [m for m in mods if m[0][1] <= cap]
    if mods:
        k, txt, how = min(mods, key=lambda m: m[0])
        # already checked by check_qasm; the CX polish of _and_finish is optional (skipped at the time limit)
        r = _stage(log, "contrib: polish", _and_finish, path, txt, k, qout + ".contrib.qasm", clean, 1.0)
        if not r:
            Path(qout + ".contrib.qasm").write_text(txt)
            r = dict(T=k[0], CX=k[1], extra_qubits=k[2], fidelity=1.0)
        if r:
            r["source"] = how; results.append((r, qout + ".contrib.qasm"))
    good = [x for x in results if x[0]["fidelity"] > 1 - 1e-9]
    if not good:
        if results:  # v0.3 behaviour: report the (failed) check
            import shutil; shutil.copy(results[0][1], qout); return results[0][0]
        return None
    r, f = min(good, key=lambda x: (x[0]["T"], x[0]["CX"], x[0]["extra_qubits"]))
    import shutil; shutil.copy(f, qout)
    return r


def _is_qc_text(text, path):
    """.qc circuit (else OpenQASM): decided by the content, by the file name only if the content does not say."""
    head = text.lstrip()
    if head.startswith("OPENQASM"):
        return False
    return head.startswith(".v") or str(path).lower().endswith(".qc")


def prove_cmd(circuit, ref, clean=False, budget=PROVE_SEC):
    """`qoptct prove`: exact proof that `circuit` (.qasm, or a .qc file) implements the input `ref` (.qc).
    -> dict(status=PROVED|REFUTED|UNKNOWN, method, detail). Needs the module qoptct.verify."""
    if not contrib_available("verify"):
        return {"status": "UNKNOWN", "method": "none", "detail": "module qoptct.verify is not installed"}
    text = Path(circuit).read_text()
    if _is_qc_text(text, circuit):
        names, hl = load_qc(circuit)
        if str(EV) not in sys.path:  # evo_lib lives in qoptct/portfolio (as in _optimize_hl)
            sys.path.insert(0, str(EV))
        import evo_lib
        text = low_to_qasm(len(names), evo_lib.lower_std(hl))
    r, err = call_contrib("verify", "prove", (text, str(Path(ref).resolve()), bool(clean), float(budget)),
                          seconds=float(budget))
    if err or not isinstance(r, dict):
        return {"status": "UNKNOWN", "method": "error", "detail": str(err or "bad result")[:300]}
    return r


_WORKDIRS = []  # (work directory of a run, keep it) -- removed by main() at the end of the run


def _umask():
    m = os.umask(0); os.umask(m); return m


def _qiskit_readable_final(fin, final, pre_cx, a, names, unitary, rep, cap):
    """The written measurement-based result has ^ / && conditions (CX post-pass, v0.3 compilers), which
    qiskit.qasm3.loads cannot read. Candidates: the rewritten result (nested single-bit ifs; proven again with
    qoptct.verify, PROVED only), the result before the CX post-pass (rewritten and proven if needed), the proven
    unitary result. The cheapest one by (T, CX) (balanced: CX <= cap) is written. -> new `final`."""
    cur = fin.read_text()
    st0 = qasm_stats(cur)
    opts = []
    for par, tag in ((True, "rewritten"), (False, "rewritten, copied corrections")):
        rw = make_qiskit_readable(cur, parity=par)
        if rw is not None and all(rw != o[1] for o in opts):
            opts.append((tag, rw, True))
    if pre_cx is not None and pre_cx != cur:
        if qiskit_readable(pre_cx):
            opts.append(("before the CX post-pass", pre_cx, False))
        else:
            for par in (True, False):
                r2 = make_qiskit_readable(pre_cx, parity=par)
                if r2 is not None and all(r2 != o[1] for o in opts):
                    opts.append(("before the CX post-pass, rewritten", r2, True))
    cand = []
    for tag, txt, need in opts:
        st = qasm_stats(txt)
        if st is None or (a.mode == "balanced" and st[2] > cap):
            continue
        cand.append(((st[1], st[2], st[0]), tag, txt, need, st))
    cand.append(((rep["t_out"], rep["cx_out"], len(names)), "unitary", None, False, None))
    info = {"from": [st0[1], st0[2]] if st0 else None}
    sec = float(PROVE_SEC if a.prove else min(PROVE_SEC, 300))
    for key, tag, txt, need, st in sorted(cand, key=lambda c: c[0][:2] + (c[4][3] if c[4] else 0,)):
        if tag == "unitary":
            write_unitary(str(fin), names, unitary)
            info.update(variant="unitary", reason="no Qiskit-readable form of the measurement-based result was proven")
            return {"format": "qc", "how": rep["method"] + " (fallback: no Qiskit-readable form of the "
                    "measurement-based result was proven)", "qiskit_rewrite": info}
        if need:
            if not contrib_available("verify"):
                continue
            pr, err = call_contrib("verify", "prove", (txt, str(Path(a.input).resolve()), bool(a.clean), sec), seconds=sec)
            status = pr.get("status") if isinstance(pr, dict) else f"error: {err}"
            info.setdefault("tried", []).append({"variant": tag, "t": st[1], "cx": st[2], "proof": status})
            if status != "PROVED":
                continue
            _PROVEN[_text_key(txt, a.input, a.clean)] = pr.get("method")
        fin.write_text(txt)
        info.update(variant=tag, to=[st[1], st[2]])
        return dict(final, how=final["how"] + f" + Qiskit-readable conditions ({tag})", qiskit_rewrite=info)
    return final  # not reached: the unitary candidate is always there


class InputError(Exception):
    """Bad command-line input (missing or unreadable file, malformed circuit): `error: ...`, exit code 3."""


EXIT_INPUT = 3  # usage / input error of `qoptct optimize` and `qoptct prove` (never 1 = REFUTED, 2 = UNKNOWN)


def _check_qc(path, what="input"):
    """Load a .qc file for the command line; any problem -> InputError with a one-line message."""
    try:
        return load_qc(path)
    except FileNotFoundError:
        raise InputError(f"{what} file not found: {path}") from None
    except IsADirectoryError:
        raise InputError(f"{what} is a directory, not a .qc file: {path}") from None
    except OSError as e:
        raise InputError(f"cannot read {what} file {path}: {e.strerror or e}") from None
    except ValueError as e:
        raise InputError(f"malformed {what} file: {e}") from None


def _check_circuit(path):
    """`qoptct prove` CIRCUIT: a .qc file, or OpenQASM 3 in the subset read by qoptct.verify."""
    try:
        text = Path(path).read_text(encoding="utf-8")
    except FileNotFoundError:
        raise InputError(f"circuit file not found: {path}") from None
    except UnicodeDecodeError:
        raise InputError(f"circuit file is not a text file: {path}") from None
    except OSError as e:
        raise InputError(f"cannot read circuit file {path}: {e.strerror or e}") from None
    if not text.strip():
        raise InputError(f"circuit file is empty: {path}")
    if _is_qc_text(text, path):
        _check_qc(path, "circuit")
        return
    try:
        from qoptct.verify._parse import parse_qasm_text
    except ImportError:
        return  # no prover installed: prove_cmd reports UNKNOWN
    try:
        parse_qasm_text(text)
    except (ValueError, KeyError, IndexError) as e:
        raise InputError(f"malformed OpenQASM circuit {path}: {e}") from None


class _ArgParser(argparse.ArgumentParser):
    def error(self, message):  # usage errors: exit code 3 (2 means UNKNOWN for `qoptct prove`)
        self.print_usage(sys.stderr)
        self.exit(EXIT_INPUT, f"{self.prog}: error: {message}\n")


def main():
    """Command line entry point (`qoptct`). Exit codes: optimize 0 = result written, 3 = usage/input error;
    prove 0 = PROVED, 1 = REFUTED, 2 = UNKNOWN, 3 = usage/input error."""
    global _END, _PROOF_RESERVE
    try:
        _main()
    except InputError as e:
        print(f"error: {e}", file=sys.stderr)
        sys.exit(EXIT_INPUT)
    finally:
        _END = None  # the overall time limit belongs to this run only (main() may be called again in-process)
        _PROOF_RESERVE = 0.0
        import shutil
        while _WORKDIRS:
            d, keep = _WORKDIRS.pop()
            if not keep:
                shutil.rmtree(d, ignore_errors=True)


def _main():
    global CONTRIB_SEC, REQUIRE_PROOF, AND_SEC
    p = _ArgParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True, parser_class=_ArgParser)
    o = sp.add_parser("optimize", help="optimize a .qc circuit (T / CX count); every written result is checked",
                      description="Optimize a Clifford+T circuit given as a Feynman .qc file. Exit code 0: result "
                                  "written (and checked); 3: usage or input error.")
    o.add_argument("input", help="input circuit (.qc)")
    o.add_argument("-o", "--output", help="output file: OpenQASM 3 if it ends in .qasm, else .qc for a unitary result (measurement-based results are always OpenQASM 3)")
    o.add_argument("--net", action="store_true"); o.add_argument("--report")
    o.add_argument("--mode", choices=["cx", "balanced", "t"], default="cx", help="cx: minimum CX, unitary (default); balanced: fewer T with CX <= --cx-factor * unitary CX; t: minimum T (measurements/ancillas allowed), higher CX")
    o.add_argument("--cx-factor", type=float, default=1.5, help="allowed CX growth in balanced mode (default 1.5)")
    o.add_argument("--clean", action="store_true", help="clean-ancilla model: .qc wires not listed in the .i line start in |0> and must end in |0>; checked with verify_qasm clean=True")
    o.add_argument("--and-variants", type=int, default=120, help="number of variants of the logical-AND compiler in balanced/t modes; 0 = off")
    o.add_argument("--time-limit", type=float, default=TIME_LIMIT, help=f"overall seconds of a run (default {TIME_LIMIT:g}, $QOPT_TIME_LIMIT; 0 = none): all searches (also the unitary portfolio) stop there and the best result checked so far is written (at worst the input); only the final check of that result (normally under a minute, at most ~15 min) may follow")
    o.add_argument("--and-sec", type=float, default=AND_SEC, help=f"seconds for the variant search of the logical-AND compiler (default {AND_SEC:g}, $QOPT_AND_SEC; 0 = unlimited)")
    o.add_argument("--sel-sec", type=float, default=600, help="seconds for the partial H-gadgetization search (balanced/t)")
    o.add_argument("--vartodd", action="store_true", help="also try an external VarTODD binary ($QOPT_VARTODD_BIN) in TODD steps (only verified results are used)")
    o.add_argument("--tgad-sec", type=float, default=120, help="seconds of CX polishing (tg3) in --tgad mode (0 = off)")
    o.add_argument("--tgad", action="store_true", help="full H gadgets + TODD with deterministic corrections; low T, much higher CX; -o is OpenQASM 3")
    o.add_argument("--post", action="store_true", help="CX post-pass (resynthesis of H-free regions; up to minutes per circuit)")
    o.add_argument("--meas", action="store_true", help="allow measurements (Gidney logical AND, clean-ancilla model of gidney.py): fewer T; -o is OpenQASM 3")
    g4 = o.add_argument_group("optional modules (v0.4)")
    for s in SOURCES:
        g4.add_argument(f"--no-{s}", action="store_true", help=f"do not use the candidate source qoptct.contrib.{s} (balanced/t)")
    g4.add_argument("--no-cxpass", action="store_true", help="do not run the CX post-pass qoptct.contrib.cxpass")
    g4.add_argument("--no-contrib", action="store_true", help="disable all optional modules (v0.3 behaviour)")
    g4.add_argument("--contrib-sec", type=float, default=CONTRIB_SEC, help=f"time limit per module call in seconds (default {CONTRIB_SEC:g}, $QOPT_CONTRIB_SEC)")
    g4.add_argument("--source-sec", action="append", default=[], metavar="NAME=SEC", help="time limit for one module (arith, bandc, gf2m, cxpass), e.g. --source-sec bandc=900; default --contrib-sec (bandc: at least 300 s above 22 qubits or 500 gates, 600 s above 200 wires or 2000 gates)")
    g4.add_argument("--budget", type=float, default=1500, help="total seconds for all optional-module stages (default 1500; 0 = unlimited)")
    g4.add_argument("--require-proof", action="store_true", help="above 22 qubits accept module results only with an exact proof (default: also when the proof is out of reach but the sampled checks of qoptct.verify pass; $QOPT_REQUIRE_PROOF=1)")
    o.add_argument("--keep-temp", action="store_true", help="keep the work directory with the intermediate files of the stages (its path is in the report as temp_dir)")
    g4.add_argument("--prove", action="store_true", help="finally prove the written result exactly with qoptct.verify (status in the report; a REFUTED result is replaced by the proven unitary one)")
    pv = sp.add_parser("prove", help="exact equivalence proof of a circuit against an input .qc (module qoptct.verify)",
                       description="Exact equivalence proof of a circuit against its input .qc. Prints JSON; exit code "
                                   "0 PROVED, 1 REFUTED, 2 UNKNOWN, 3 usage or input error.")
    pv.add_argument("circuit", help="optimized circuit: OpenQASM 3 (as written by qoptct) or .qc")
    pv.add_argument("reference", help="input circuit (.qc)")
    pv.add_argument("--clean", action="store_true", help="clean-ancilla model (wires not in the .i line start and end in |0>)")
    pv.add_argument("--budget", type=float, default=PROVE_SEC, help=f"seconds (default {PROVE_SEC})")
    sp.add_parser("kprove", add_help=False,
                  help="exact BDD-based proof for large classical-reference circuits (qoptct.verify.kprove); "
                       "see `qoptct kprove --help`")
    sp.add_parser("version")
    if sys.argv[1:2] == ["kprove"]:
        from .verify.kprove import main as kprove_main
        sys.exit(kprove_main(sys.argv[2:]))
    a = p.parse_args()
    if a.cmd == "version":
        print(f"qoptct {__version__}"); return
    if a.cmd == "prove":
        _check_qc(a.reference, "reference")
        _check_circuit(a.circuit)
        r = prove_cmd(a.circuit, a.reference, a.clean, a.budget)
        print(json.dumps(r, ensure_ascii=False))
        det = str(r.get("detail", ""))
        if r.get("status") not in ("PROVED", "REFUTED") and (r.get("method") == "parse" or "unsupported gate" in det):
            raise InputError(f"the prover cannot read the input: {det}"[:300])
        sys.exit({"PROVED": 0, "REFUTED": 1}.get(r.get("status"), 2))
    _check_qc(a.input, "input")
    if a.output:
        od = Path(a.output).resolve().parent
        if not od.is_dir():
            raise InputError(f"output directory does not exist: {od}")
        if Path(a.output).is_dir():
            raise InputError(f"output is a directory: {a.output}")
    if a.net:
        try:
            require_torch()
        except ImportError as e:
            p.exit(2, f"{p.prog}: error: {e}\n")  # optional dependency missing (unchanged since 0.3)
    CONTRIB_SEC = a.contrib_sec
    AND_SEC = a.and_sec
    for kv in a.source_sec:
        try:
            k, v = kv.split("="); MODULE_SEC[k.strip()] = float(v)
        except ValueError:
            p.error(f"--source-sec expects NAME=SEC, got {kv!r}")
    REQUIRE_PROOF = REQUIRE_PROOF or a.require_proof
    use_contrib = not a.no_contrib
    sources = tuple(s for s in SOURCES if use_contrib and not getattr(a, f"no_{s}") and contrib_available(f"contrib.{s}"))
    use_cxpass = use_contrib and not a.no_cxpass and contrib_available("contrib.cxpass")
    budget, mlog = Budget(a.budget), []
    import tempfile, shutil
    tmpdir = None
    if not a.output:  # the result is always written (to a temporary file without -o) so the report describes it
        tmpdir = tempfile.mkdtemp(prefix="qopt_out_")
        a.output = str(Path(tmpdir) / (Path(a.input).stem + ".out"))
    # all intermediate files (stage results, *.and.qasm, *.sel.qasm, *_tg3.qasm, ...) live in a private work
    # directory; the final result is moved to -o atomically at the end (--keep-temp keeps the work directory)
    final_out = None if tmpdir else a.output
    work = tempfile.mkdtemp(prefix="qopt_work_")
    _WORKDIRS.append((work, a.keep_temp))
    a.output = str(Path(work) / ("result" + Path(a.output).suffix))
    t = time.time()
    global _END
    _END = t + a.time_limit if a.time_limit and a.time_limit > 0 else None
    global _PROOF_RESERVE
    _PROOF_RESERVE = min(float(PROVE_SEC), 0.25 * a.time_limit) if _END is not None else 0.0
    names, ref, best, how = optimize(a.input, a.net)
    budget = Budget(a.budget)  # --budget covers the optional-module stages only: its clock starts after the portfolio
    import evo_lib
    if a.post and how:
        sys.path.insert(0, str(ROOT / "contrib" / "eng2"))
        spec = importlib.util.spec_from_file_location("eng2_improve", ROOT / "contrib" / "eng2" / "improve.py")
        e2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(e2)
        b2 = _stage(mlog, "post", e2.improve, names, best, ref, default=best)  # its input if no improvement is proven
        if evo_lib.key(b2) < evo_lib.key(best) and _stage(mlog, "post: proof", prove, b2, ref, len(names)) is True:
            best, how = b2, how + "+post"
    best0, cx0 = best, evo_lib.cx(best)  # portfolio result: start of the H-gadget searches and base of the balanced cap
    if use_cxpass and how:
        c = cxpass_unitary(best, ref, len(names), budget, mlog, path=a.input)
        if c is not None:
            best, how = c, how + "+cxpass"
    rep = {"input": a.input, "qubits": len(names), "cx_in": evo_lib.cx(ref), "t_in": evo_lib.tcount(ref),
           "cx_out": evo_lib.cx(best), "t_out": evo_lib.tcount(best), "method": how or "input (not proven)",
           "verified": how is not None, "seconds": round(time.time() - t, 1)}
    segs = None
    final = {"format": "qc", "how": rep["method"]}
    if a.meas:
        _, segs, T, CX, ok, mhow = _stage(mlog, "meas", optimize_meas, a.input,
                                          default=(None, None, None, None, None, "stopped at the time limit"))
        rep["meas"] = {"t": T, "cx": CX, "verified": ok, "proof": mhow}
        if not (ok is True and T is not None and T < rep["t_out"]):
            segs = None  # not better without measurements or not proven: keep the unitary result
        rep["meas"]["used"] = segs is not None
        if segs is not None:
            final = {"format": "qasm", "how": "meas: " + str(mhow)}
    if a.mode == "t":
        a.tgad = True
    if a.vartodd:
        sys.path.insert(0, str(ROOT / "contrib" / "eng1"))
        import vartodd_bridge
        rep["vartodd"] = vartodd_bridge.install()
    out0 = a.output
    # the v0.3 measurement-based stages (partial/full H gadgets, logical-AND compilers) accept only results that the
    # simulator can check (<= SIM_MAX_QUBITS qubits in total); on larger inputs they cannot succeed and are skipped
    small = len(names) <= SIM_MAX_QUBITS
    # v0.4 candidate sources first (bounded by --budget): the v0.3 stages below are open-ended searches that the
    # overall time limit may cut; the final choice is the minimum over all candidates, so the order does not matter
    pre_mods = contrib_candidates(a.input, a.clean, sources, budget, mlog) \
        if a.mode in ("balanced", "t") and sources else []
    if a.mode in ("balanced", "t") and not small:
        rep["skipped_v03_stages"] = f"input has {len(names)} > {SIM_MAX_QUBITS} qubits: H-gadget and v0.3 AND stages skipped"
    if a.mode == "balanced" and small:
        qf = str(Path(a.output) if a.output else Path("/tmp") / (Path(a.input).stem + "_bal.qasm"))
        try:
            r = _stage(mlog, "balanced: H gadgets", run_sel, a.input, names, best0, a.cx_factor * cx0, a.sel_sec, qf)
        except Exception as e:
            r = None; rep["balanced_error"] = str(e)[:200]
        ok = bool(r and r.get("fidelity") and r["fidelity"] > 1 - 1e-9 and r["T"] < rep["t_out"])
        rep["balanced"] = {"t": r and r.get("T"), "cx": r and r.get("CX"), "extra_qubits": r and r.get("extra_qubits"),
                           "independent_fidelity": r and r.get("fidelity"), "used": ok and bool(a.output), "qasm": qf}
        if ok and a.output:
            a.output = None; final = {"format": "qasm", "how": f"partial H gadgets, simulator F={r['fidelity']:.12f}"}
    if a.tgad and small:
        spec = importlib.util.spec_from_file_location("eng1_det", ROOT / "contrib" / "eng1" / "det.py")
        sys.path.insert(0, str(ROOT / "contrib" / "eng1"))
        d1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d1)
        qf = str(Path(a.output) if a.output else Path("/tmp") / (Path(a.input).stem + "_tgad.qasm"))
        r = _stage(mlog, "t: H gadgets + TODD", lambda: d1.improve(a.input, qasm_out=qf, mode="B")) or {}
        import verify_qasm
        fid = _stage(mlog, "t: H gadgets + TODD check", lambda: limited(PROVE_SEC, verify_qasm.check, qf, a.input)) \
            if r and r.get("extra_qubits", 99) + len(names) <= 22 else None
        rep["tgad"] = {"t": r.get("T"), "cx": r.get("CX"), "extra_qubits": r.get("extra_qubits"), "measurements": r.get("measurements"),
                       "agent_verified": r.get("verified"), "independent_fidelity": fid, "qasm": qf}
        # CX polishing by agent tg3 (T unchanged); used only after the independent check
        if fid is not None and fid > 1 - 1e-9 and a.tgad_sec > 0:
            try:
                sys.path.insert(0, str(ROOT / "contrib" / "tg3"))
                spec3 = importlib.util.spec_from_file_location("tg3_opt", ROOT / "contrib" / "tg3" / "tg3_opt.py")
                g3 = importlib.util.module_from_spec(spec3); spec3.loader.exec_module(g3)
                q3 = qf[:-5] + "_tg3.qasm" if qf.endswith(".qasm") else qf + "_tg3"
                r3 = _stage(mlog, "t: tg3 polish", lambda: g3.improve(a.input, qasm_out=q3, seconds=a.tgad_sec, det_qasm=qf,
                                                                       verify=False, log=lambda *x: None))
                f3 = _stage(mlog, "t: tg3 check", lambda: limited(PROVE_SEC, verify_qasm.check, q3, a.input)) if r3 else None
                if r3 is None or f3 is None:
                    raise RuntimeError("tg3: no checked result (time limit)")
                rep["tgad"]["tg3"] = {"t": r3.get("T"), "cx": r3.get("CX"), "independent_fidelity": f3}
                if f3 > 1 - 1e-9 and r3.get("T") == r.get("T") and r3.get("CX", 1e9) < (r.get("CX") or 1e9):
                    shutil.copy(q3, qf)
                    rep["tgad"].update({"cx": r3.get("CX"), "independent_fidelity": f3, "post": "tg3"})
            except (_TO, Exception) as e:
                rep["tgad"]["tg3_error"] = str(e)[:200]
        if a.output and not (fid is not None and fid > 1 - 1e-9):
            write_unitary(a.output, names, best); rep["tgad"]["used"] = False
        else:
            rep["tgad"]["used"] = bool(a.output)
            if a.output:
                final = {"format": "qasm", "how": f"H gadgets + TODD, simulator F={rep['tgad']['independent_fidelity']:.12f}"}
        a.output = None if rep["tgad"].get("used") else a.output
        if a.mode == "t":
            qs = str(Path(out0 or "/tmp/" + Path(a.input).stem + "_t.qasm")) + ".sel.qasm"
            try:
                r2 = _stage(mlog, "t: H gadgets", run_sel, a.input, names, best0, 1e9, a.sel_sec, qs)
            except Exception as e:
                r2 = None; rep["tgad"]["sel_error"] = str(e)[:200]
            if r2:
                rep["tgad"]["sel"] = {"t": r2.get("T"), "cx": r2.get("CX"), "extra_qubits": r2.get("extra_qubits"), "independent_fidelity": r2.get("fidelity")}
                cur = (rep["tgad"]["t"], rep["tgad"]["cx"]) if rep["tgad"].get("used") or (not out0 and rep["tgad"].get("independent_fidelity", 0) and rep["tgad"]["independent_fidelity"] > 1 - 1e-9) else (rep["t_out"], 10**9)
                if r2["fidelity"] and r2["fidelity"] > 1 - 1e-9 and (r2["T"], r2["CX"]) < tuple(x if x is not None else 10**9 for x in cur):
                    rep["tgad"]["choice"] = "sel"
                    if out0:
                        shutil.copy(qs, out0); a.output = None; rep["tgad"]["used"] = True
                        final = {"format": "qasm", "how": f"partial H gadgets, simulator F={r2['fidelity']:.12f}"}
    cap = 1e9
    if a.mode in ("balanced", "t") and (a.and_variants > 0 or sources):
        if a.mode == "balanced":
            b = rep.get("balanced", {})
            chosen = (b["t"], b["cx"]) if b.get("t") and b.get("independent_fidelity") and b["independent_fidelity"] > 1 - 1e-9 and b["t"] < rep["t_out"] else (rep["t_out"], rep["cx_out"])
            cap = a.cx_factor * cx0
        else:
            g = rep.get("tgad", {})
            if g.get("choice") == "sel":
                chosen = (g["sel"]["t"], g["sel"]["cx"])
            elif g.get("independent_fidelity") and g["independent_fidelity"] > 1 - 1e-9:
                chosen = (g["t"], g["cx"])
            else:
                chosen = (rep["t_out"], rep["cx_out"])
        qa = str(Path(out0 or "/tmp/" + Path(a.input).stem + "_and.qasm")) + ".and.qasm"
        try:
            ra = run_and(a.input, names, cap, a.and_variants if small else 0, qa, clean=a.clean, sources=sources, budget=budget,
                         log=mlog, mods=pre_mods)
        except Exception as e:
            ra = None; rep["and_error"] = str(e)[:200]
        if ra:
            use = ra["fidelity"] > 1 - 1e-9 and (ra["T"], ra["CX"]) < chosen
            rep["and"] = {"t": ra["T"], "cx": ra["CX"], "extra_qubits": ra["extra_qubits"], "independent_fidelity": ra["fidelity"],
                          "chosen": use, "model": "clean" if a.clean else "strict", "source": ra.get("source")}
            if use and out0:
                shutil.copy(qa, out0); a.output = None
                final = {"format": "qasm", "how": f"{ra.get('source')}" + (f" + {ra['post']}" if ra.get("post") else "")}
    if a.output:
        if segs is not None:
            save_qasm3(a.output, len(names), segs)
        else:
            write_unitary(a.output, names, best)  # OpenQASM 3 for -o *.qasm, else .qc
    # ---- v0.4: CX post-pass of a measurement-based result, final exact proof, description of the written result
    fin = Path(out0)
    pre_cx = fin.read_text() if final["format"] == "qasm" else None  # the result before the CX post-pass
    if final["format"] == "qasm" and use_cxpass:
        if a.mode == "balanced":
            better = lambda s0, s1: (s1[1], s1[2]) < (s0[1], s0[2]) and s1[2] <= cap
        else:
            better = lambda s0, s1: (s1[1], s1[2]) < (s0[1], s0[2])
        txt = cxpass_qasm(fin.read_text(), a.input, a.clean, budget, mlog, better)
        if txt is not None:
            fin.write_text(txt); final["how"] += " + cxpass"
    # the searches end here; the final check of the written result is not cut by --time-limit (see --help)
    rep["search_seconds"] = round(time.time() - t, 1)
    _END = None
    if final["format"] == "qasm" and not qiskit_readable(fin.read_text()):
        final = _qiskit_readable_final(fin, final, pre_cx, a, names, best if how else ref, rep, cap)
        if final["format"] == "qc":
            best = best if how else ref
    if a.prove or (final["format"] == "qasm" and contrib_available("verify")):
        # measurement-based results are always re-checked exactly when qoptct.verify is installed (short budget;
        # --prove: full budget, also for unitary results); a REFUTED result is never written
        want = PROVE_SEC if a.prove else SMALL_PROVE_SEC
        sec = want if budget.end is None else max(0.0, min(want, budget.end - time.time()))  # not cut by --time-limit
        pr = prove_cmd(str(fin), a.input, a.clean, max(sec, 30 if a.prove else 10))
        final["proof"] = {k: pr[k] for k in ("status", "method", "detail", "time", "counterexample") if k in pr}
        earlier = _PROVEN.get(_text_key(fin.read_text(), a.input, a.clean)) if final["format"] == "qasm" else None
        if pr.get("status") == "UNKNOWN" and earlier:
            # the final proof ran out of its short budget, but this exact text was proven (PROVED) earlier in the run
            final["check"] = f"PROVED earlier in this run ({earlier}); final re-proof: UNKNOWN (time)"
        elif final["format"] == "qasm" and pr.get("status") == "UNKNOWN":
            # no exact proof: the written measurement-based result must pass the strong sampled check
            ok, h2 = strong_check(fin.read_text(), a.input, a.clean)
            final["check"] = h2
            if not ok:
                pr = dict(pr, status="REFUTED", detail=f"strong sampled check failed: {h2}", not_confirmed=True)
        elif pr.get("status") == "PROVED":
            final["check"] = f"PROVED ({pr.get('method')})"
        if pr.get("status") == "REFUTED":  # never hand out a refuted circuit: fall back to the proven unitary one
            write_unitary(str(fin), names, best if (how and final["format"] == "qasm") else ref)
            best = best if (how and final["format"] == "qasm") else ref; how = how if final["format"] == "qasm" else None
            why = "the final check could not confirm the result" if pr.get("not_confirmed") else "the result was REFUTED"
            final = {"format": "qc", "how": (rep["method"] if how else "input") + f" (fallback: {why})", "proof": pr}
    if final["format"] == "qasm":
        txt = fin.read_text()
        st = qasm_stats(txt) or (int(re.search(r"qubit\[(\d+)\]", txt).group(1)), len(re.findall(r"\b(t|tdg)\s+q", txt)),
                                 len(re.findall(r"\b(cx|cz)\s+q", txt)), len(re.findall(r"measure", txt)))
        final.update(t=st[1], cx=st[2], extra_qubits=st[0] - len(names), measurements=st[3], verified=True)
    else:
        final.update(t=evo_lib.tcount(best), cx=evo_lib.cx(best), extra_qubits=0, measurements=0, verified=how is not None)
        if how is None:
            final["note"] = "no improvement was proven (or the time limit stopped the search): the input is written"
    final["file"] = "qc" if fin.read_text()[:400].lstrip().startswith(".v") else "openqasm3"
    if final_out and final["file"] == "openqasm3" and not str(final_out).lower().endswith(QASM_SUFFIXES):
        print(f"warning: {final_out}: the result has measurements and is written as OpenQASM 3 (use a .qasm file "
              "name)", file=sys.stderr)
    if final_out:  # atomic: a temporary file next to -o, then rename (never a half-written or intermediate -o)
        dst = Path(final_out)
        fd, tmpf = tempfile.mkstemp(prefix="." + dst.name + ".", suffix=".tmp", dir=str(dst.resolve().parent))
        os.close(fd)
        shutil.copyfile(fin, tmpf)
        os.chmod(tmpf, 0o644 & ~_umask())
        os.replace(tmpf, dst)
    if a.keep_temp:
        rep["temp_dir"] = work
    rep["final"] = final
    rep["modules"] = {"sources": list(sources), "cxpass": use_cxpass, "verify": contrib_available("verify"), "log": mlog}
    rep["seconds_total"] = round(time.time() - t, 1)
    if tmpdir:
        shutil.rmtree(tmpdir, ignore_errors=True)
    print(json.dumps(rep, ensure_ascii=False))
    if a.report:
        Path(a.report).write_text(json.dumps(rep, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
