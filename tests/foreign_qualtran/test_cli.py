"""``qoptct qualtran``: exit codes and output.  Exit code 0 only for PROVED."""
import json
import subprocess
import sys

import pytest

from qoptct.foreign.qualtran import cli


def run(args, capsys):
    code = cli.main(args)
    out = capsys.readouterr()
    return code, out.out, out.err


def test_results_table(capsys):
    code, out, _ = run(["results"], capsys)
    assert code == 0
    for name in ("montgomery_mul", "kaliski_inverse", "pair_montgomery", "pair_kaliski", "ecadd", "ecwindowaddr"):
        assert name in out
    assert "PROVED" in out and "PROVED_OUTSIDE_G1_G3" in out


def test_results_filter_json_long(capsys):
    code, out, _ = run(["results", "--name", "ecwindowaddr", "--n", "256", "--json"], capsys)
    es = json.loads(out)
    assert code == 0 and {e["w"] for e in es} == {1, 2, 4} and all(e["name"] == "ecwindowaddr" for e in es)
    code, out, _ = run(["results", "--name", "kaliski_inverse", "--n", "16", "--long"], capsys)
    assert code == 0 and "reproduce:  qoptct qualtran prove kaliski --n 16 --p 32749" in out and "assumptions:" in out
    code, _, err = run(["results", "--name", "ecadd", "--n", "17"], capsys)
    assert code == cli.EXIT_USAGE and "no such result" in err


def test_usage_errors(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["prove", "nonsense", "--n", "4"])
    assert e.value.code == 2                                      # argparse
    capsys.readouterr()
    code, _, err = run(["prove", "ecadd", "--n", "20"], capsys)  # no default modulus at n = 20 (or no Qualtran: also exit 3)
    assert code == cli.EXIT_USAGE and err.startswith("error:")
    code, _, err = run(["prove", "montgomery", "--n", "16", "--p", "__import__('os')"], capsys)
    assert code == cli.EXIT_USAGE


def test_modulus_parser():
    assert cli._int_expr("2**255-19") == 2 ** 255 - 19 and cli._int_expr("32749") == 32749
    for bad in ("os.system('x')", "1e3", "", "9**9**9**9**9" + "9" * 80):
        with pytest.raises(cli.UsageError):
            cli._int_expr(bad)


def test_module_entry_point_and_main_cli():
    # `python -m` of the module and the `qoptct qualtran` subcommand reach the same code
    r = subprocess.run([sys.executable, "-B", "-m", "qoptct.foreign.qualtran.cli", "results", "--name", "pair_kaliski"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and "pair_kaliski" in r.stdout
    r = subprocess.run([sys.executable, "-B", "-c", "import sys; sys.argv=['qoptct','qualtran','results','--name','ecadd','--n','16']; "
                        "from qoptct.qopt_core import main; main()"], capture_output=True, text=True)
    assert r.returncode == 0 and "ecadd" in r.stdout, r.stderr


# --- the commands that need Qualtran (and Z3 for the proofs)
def _need():
    pytest.importorskip("qualtran")
    pytest.importorskip("z3")


def test_prove_montgomery_n16_w1(capsys):
    _need()
    code, out, _ = run(["prove", "montgomery", "--n", "16", "--w", "1"], capsys)
    assert code == 0 and out.startswith("PROVED") and "program hash" in out


def test_check_command(capsys):
    pytest.importorskip("qualtran")
    code, out, _ = run(["check", "qualtran.bloqs.arithmetic.addition.Add(QUInt(3))"], capsys)
    assert code == 0 and out.startswith("PROVED"), out
    code, _, err = run(["check", "no_such_name(3)"], capsys)
    assert code == cli.EXIT_USAGE and "cannot build the bloq" in err


def test_original_ecadd_exit_code_is_not_zero(capsys):
    """The original scheme is proved only outside the groups G1-G3: that is a different statement, the exit code is 1."""
    _need()
    code, out, _ = run(["prove", "ecadd", "--n", "16", "--original"], capsys)
    assert code == cli.EXIT_NOT_PROVED and out.startswith("PROVED_OUTSIDE_G1_G3")
