"""Throwaway runner that mirrors examples/repair_cascade.ipynb.

Run from the repo root so the local packages are importable:

    virtenv/bin/python examples/_run_repair_cascade.py
"""

import contextlib
import io
import os
import sys
import time

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from systems.lindivs import LinDivs, OrdLinDivs
from tests.test_norm_xval import (
    _with_timeout, _NormTimeout, encode_leaf, _check,
)


def repair_cascade(n: int) -> LinDivs:
    assert n >= 2
    veclen = (n + 1) + 1

    def vec(pairs, const=0):
        v = [0] * veclen
        for i, c in pairs:
            v[i] = c
        v[-1] = const
        return tuple(v)

    F, G = [], []
    F.append(vec([(n, 1)]))
    G.append(vec([(n - 1, 3)]))
    for k in range(2, n + 1):
        F.append(vec([(k - 1, 1), (k - 2, 1)]))
        G.append(vec([(n, 1)]))
    return LinDivs(tuple(F), tuple(G))


def lipshitz_coef_metrics(leaf):
    F, G = leaf.F, leaf.G
    coeffs = [c for row in F for c in row] + [c for row in G for c in row]
    # Tautological residual: every row of F equals the matching row of G,
    # so each constraint reduces to "f | f" (trivially satisfied).
    tautological = (len(F) > 0 and tuple(F) == tuple(G))
    if not coeffs:
        return {"n_F": 0, "n_G": 0, "max_abs": 0,
                "max_bits": 0, "total_bits": 0, "tautological": False}
    max_abs = int(max(abs(c) for c in coeffs))
    return {
        "n_F": len(F),
        "n_G": len(G),
        "max_abs": max_abs,
        "max_bits": max_abs.bit_length(),
        "total_bits": sum(int(abs(c)).bit_length() for c in coeffs),
        "tautological": tautological,
    }


def P(i: int) -> int:
    p = 3
    for _ in range(i):
        p = p * (p - 1)
    return p


def analyse(n, timeout_s=120, check_sym_inc=True, use_all_cx_inc=False):
    lds = repair_cascade(n)
    t0 = time.time()
    try:
        sink = io.StringIO()
        with contextlib.redirect_stdout(sink):
            leaves = _with_timeout(timeout_s, lds.norm,
                                   check_sym_inc=check_sym_inc,
                                   use_all_cx_inc=use_all_cx_inc)
    except _NormTimeout:
        return {"n": n, "status": "TIMED_OUT",
                "timeout_s": timeout_s, "wall_s": time.time() - t0}
    except AssertionError as e:
        import traceback as _tb
        return {"n": n, "status": "NORM_BUG",
                "error": "AssertionError in LinDivs.norm()",
                "traceback": _tb.format_exc(),
                "wall_s": time.time() - t0}
    except Exception as e:
        return {"n": n, "status": "ERROR",
                "error": repr(e), "wall_s": time.time() - t0}
    wall = time.time() - t0

    n_sat = n_unsat = n_unknown = 0
    sat_trivial = 0       # SAT, fully substituted (n_F == 0)
    sat_taut = 0          # SAT, residual but every row tautological (F==G)
    sat_nontriv = []      # (leaf, metrics): residual AND non-tautological
    for leaf in leaves:
        r = _check(encode_leaf(leaf))
        if r is True:
            n_sat += 1
            if isinstance(leaf, OrdLinDivs):
                m = lipshitz_coef_metrics(leaf)
            else:
                m = {"n_F": 0, "n_G": 0, "max_abs": 0,
                     "max_bits": 0, "total_bits": 0, "tautological": False}
            if m["n_F"] == 0:
                sat_trivial += 1
            elif m["tautological"]:
                sat_taut += 1
            else:
                sat_nontriv.append((leaf, m))
        elif r is False:
            n_unsat += 1
        else:
            n_unknown += 1

    nt_min = (min(sat_nontriv, key=lambda lm: lm[1]["max_abs"])
              if sat_nontriv else None)
    nt_max = (max(sat_nontriv, key=lambda lm: lm[1]["max_abs"])
              if sat_nontriv else None)

    return {
        "n": n,
        "status": "ok",
        "settings": {"check_sym_inc": check_sym_inc,
                     "use_all_cx_inc": use_all_cx_inc},
        "n_leaves": len(leaves),
        "n_sat": n_sat, "n_unsat": n_unsat, "n_unknown": n_unknown,
        "n_sat_trivial": sat_trivial,
        "n_sat_tautological": sat_taut,
        "n_sat_nontrivial": len(sat_nontriv),
        "P_floor": P(n - 1),
        "P_cap": P(n),
        "nt_min_max_abs": nt_min[1]["max_abs"] if nt_min else None,
        "nt_max_max_abs": nt_max[1]["max_abs"] if nt_max else None,
        "nt_min_max_bits": nt_min[1]["max_bits"] if nt_min else None,
        "nt_max_max_bits": nt_max[1]["max_bits"] if nt_max else None,
        "nt_min_leaf": nt_min[0] if nt_min else None,
        "nt_max_leaf": nt_max[0] if nt_max else None,
        "wall_s": wall,
    }


def report(rows):
    print()
    print("n | leaves | SAT | empty | taut | n-triv | nt min |c| | nt max |c| | P_{n-1}    | verdict")
    print("--+--------+-----+-------+------+--------+------------+------------+------------+--------")
    compact_witnesses = []
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['n']} | {r['status']}  (wall {r['wall_s']:.1f}s)")
            continue
        n_nt = r["n_sat_nontrivial"]
        nmin = r["nt_min_max_abs"]
        nmax = r["nt_max_max_abs"]
        p_floor = r["P_floor"]
        if r["n_sat"] == 0:
            verdict = "UNSAT (no SAT branches)"
        elif n_nt == 0:
            verdict = "all SAT branches collapse (empty or tautological residual)"
        elif nmin is not None and nmin < p_floor:
            verdict = f"COMPACT NT max|c|={nmin} < P_{r['n']-1}={p_floor}"
            compact_witnesses.append(r)
        else:
            verdict = f"all NT SAT >= P_{r['n']-1}={p_floor}"
        print(
            f"{r['n']} | {r['n_leaves']:>6} | {r['n_sat']:>3} | "
            f"{r['n_sat_trivial']:>5} | {r['n_sat_tautological']:>4} | "
            f"{n_nt:>6} | {str(nmin):>10} | {str(nmax):>10} | "
            f"{p_floor:>10} | {verdict}"
        )

    print()
    if compact_witnesses:
        print("=== COMPACT NON-TRIVIAL RESIDUAL WITNESSES ===")
        for r in compact_witnesses:
            print(f"\nn = {r['n']}, non-trivial residual max|c| = {r['nt_min_max_abs']}:")
            print(r["nt_min_leaf"])
    else:
        print("=== NO COMPACT NON-TRIVIAL RESIDUAL ESCAPE ===")
        print("Every conclusive SAT branch either (a) substituted to empty "
              "residual, (b) reduced to tautologies (F == G), or "
              "(c) carried at least one coefficient meeting the cascade floor.")


def main():
    TIMEOUTS = {2: 60, 3: 60, 4: 120, 5: 240}
    rows = []
    for n in (2, 3, 4, 5):
        print(f"--- n = {n} (timeout {TIMEOUTS[n]}s) ---", flush=True)
        r = analyse(n, timeout_s=TIMEOUTS[n])
        if r["status"] == "ok":
            print(
                f"  leaves={r['n_leaves']} (sat={r['n_sat']} "
                f"[empty={r['n_sat_trivial']}, taut={r['n_sat_tautological']}, "
                f"nontriv={r['n_sat_nontrivial']}], "
                f"unsat={r['n_unsat']}, unknown={r['n_unknown']})\n"
                f"  P_floor=P_{n - 1}={r['P_floor']}  P_cap=P_{n}={r['P_cap']}\n"
                f"  non-trivial residual max_abs: "
                f"min={r['nt_min_max_abs']}, max={r['nt_max_max_abs']} "
                f"(bits {r['nt_min_max_bits']}..{r['nt_max_max_bits']})\n"
                f"  wall={r['wall_s']:.2f}s",
                flush=True,
            )
        else:
            print(f"  {r['status']} after {r['wall_s']:.1f}s", flush=True)
            if r["status"] == "NORM_BUG":
                print("  (LinDivs.norm() AssertionError -- "
                      "pre-existing bug in normalizer)", flush=True)
        rows.append(r)
    report(rows)


if __name__ == "__main__":
    main()
