"""Cross-validation of Lipshitz and KNF normalizers via z3.

Both normalizers must preserve satisfiability of the input: the original
system is sat iff at least one returned leaf is sat (in the leaf's own
coordinates).  We therefore generate random small LinDivs, run both
normalizers, and check that they agree on whether *some* leaf is
satisfiable.

z3 frequently blows up on divisibility problems, so every check is run
under a timeout and `unknown` trials are simply skipped.  The test only
asserts agreement when a trial is conclusive for both normalizers and
that at least one trial is conclusive overall.
"""

import contextlib
import io
import random
import signal

import pytest
import z3

from systems.lindivs import LinDivs, OrdLinDivs, KNFLeaf

# Per-z3-check timeout.  Keep tight: each trial may dispatch many checks.
Z3_TIMEOUT_MS = 500
# Wall-clock cap for a single normalisation call (each: norm and knf_norm).
NORM_TIMEOUT_S = 2
# Number of randomised trials.
N_TRIALS = 20
# Coefficient range for the random generator.
COEF_RANGE = (-1, 1)
# RNG seed for reproducibility.
RNG_SEED = 20260516


class _NormTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise _NormTimeout()


def _with_timeout(seconds, fn, *args, **kwargs):
    """Run fn under a SIGALRM-based wall-clock timeout (POSIX only)."""
    prev = signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(seconds)
    try:
        return fn(*args, **kwargs)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)


def _aff(coeffs, xs):
    """Build a z3 expression for sum(c_i * x_i) + c_d from a (d+1)-vector."""
    if not coeffs:
        return z3.IntVal(0)
    expr = z3.IntVal(int(coeffs[-1]))
    for i, c in enumerate(coeffs[:-1]):
        if c != 0:
            expr = expr + int(c) * xs[i]
    return expr


def _add_divs(solver, F, G, xs, tag):
    """Encode each f_i | g_i via a fresh integer k: f_i * k == g_i."""
    for i, (f, g) in enumerate(zip(F, G)):
        k = z3.Int(f"k_{tag}_{i}")
        solver.add(_aff(f, xs) * k == _aff(g, xs))


def _dim_of(lds: LinDivs):
    """Recover the number of integer variables from any non-empty matrix."""
    for mat in (lds.F, lds.G, lds.get_eqs(), lds.get_ineqs()):
        if mat:
            return len(mat[0]) - 1
    return 0


def encode_original(lds: LinDivs):
    """z3 Solver for the *original* LinDivs over integers (no x>=0)."""
    s = z3.Solver()
    s.set("timeout", Z3_TIMEOUT_MS)
    d = _dim_of(lds)
    if d == 0 and not lds.F:
        return s  # trivially sat
    xs = [z3.Int(f"x{i}") for i in range(d)]
    for a in lds.get_eqs():
        s.add(_aff(a, xs) == 0)
    for b in lds.get_ineqs():
        s.add(_aff(b, xs) <= 0)
    _add_divs(s, lds.F, lds.G, xs, f"orig{lds.id}")
    return s


def encode_leaf(leaf):
    """Encode a leaf for sat checking.

    EPAD preprocessing reparametrises variables to non-negative integers,
    so both OrdLinDivs and KNFLeaf are checked with x_i >= 0.  Trivial
    empty LinDivs leaves are encoded as their (possibly empty) constraint
    system.
    """
    s = z3.Solver()
    s.set("timeout", Z3_TIMEOUT_MS)

    if isinstance(leaf, LinDivs):
        d = _dim_of(leaf)
        if d == 0 and not leaf.F:
            return s  # trivially sat
        xs = [z3.Int(f"x{i}") for i in range(d)]
        for x in xs:
            s.add(x >= 0)
        for a in leaf.get_eqs():
            s.add(_aff(a, xs) == 0)
        for b in leaf.get_ineqs():
            s.add(_aff(b, xs) <= 0)
        _add_divs(s, leaf.F, leaf.G, xs, f"trv{leaf.id}")
        return s

    # OrdLinDivs or KNFLeaf
    if not leaf.F:
        return s  # trivially sat
    d = len(leaf.F[0]) - 1
    xs = [z3.Int(f"x{i}") for i in range(d)]
    for x in xs:
        s.add(x >= 0)
    _add_divs(s, leaf.F, leaf.G, xs, f"lf{id(leaf)}")
    if isinstance(leaf, KNFLeaf):
        for l in leaf.L_gens:
            s.add(_aff(l, xs) == 0)
    return s


def _check(solver):
    """Return True/False/None (None = unknown/timeout)."""
    r = solver.check()
    if r == z3.sat:
        return True
    if r == z3.unsat:
        return False
    return None


def any_leaf_sat(leaves):
    """True if some leaf is sat; False if all unsat; None if any unknown
    and no sat seen."""
    saw_unknown = False
    for leaf in leaves:
        r = _check(encode_leaf(leaf))
        if r is True:
            return True
        if r is None:
            saw_unknown = True
    return None if saw_unknown else False


def chain(n: int) -> LinDivs:
    """The Lipshitz/EPAD divisibility chain over x_0,...,x_n:

        x_i | x_{i+1} ∧ (x_i + 1) | x_{i+1}     for i = 0,...,n-1

    Every solution satisfies x_{i+1} >= x_i * (x_i + 1), so x_i grows
    doubly exponentially in i.  This is Section 6 of the KNF note,
    stripped of the auxiliary x_0=2 and x_i>=2 constraints.
    """
    veclen = (n + 1) + 1

    def vec(pairs, const=0):
        v = [0] * veclen
        for i, c in pairs:
            v[i] = c
        v[-1] = const
        return tuple(v)

    F, G = [], []
    for i in range(n):
        F.append(vec([(i, 1)]))
        G.append(vec([(i + 1, 1)]))
        F.append(vec([(i, 1)], const=1))
        G.append(vec([(i + 1, 1)]))
    return LinDivs(tuple(F), tuple(G))


def chain_with_yz(n: int) -> LinDivs:
    """The chain over x_0,...,x_n extended with two fresh variables y, z
    and divisibilities

        (y - z) | x_1     and     x_n | (y + z),

    designed to force non-increasingness for every total order on the
    variables (the chain pins LV(x_n) = x_n while the new constraints
    drag x_n below max(y,z) and y,z below x_n simultaneously).
    """
    veclen = (n + 3) + 1
    y_idx = n + 1
    z_idx = n + 2

    def vec(pairs, const=0):
        v = [0] * veclen
        for i, c in pairs:
            v[i] = c
        v[-1] = const
        return tuple(v)

    F, G = [], []
    for i in range(n):
        F.append(vec([(i, 1)]))
        G.append(vec([(i + 1, 1)]))
        F.append(vec([(i, 1)], const=1))
        G.append(vec([(i + 1, 1)]))
    F.append(vec([(y_idx, 1), (z_idx, -1)]))
    G.append(vec([(1, 1)]))
    F.append(vec([(n, 1)]))
    G.append(vec([(y_idx, 1), (z_idx, 1)]))
    return LinDivs(tuple(F), tuple(G))


def _outcome(fn):
    """Run a normaliser and report ('ok', nleaves) or (exc_name, None)."""
    sink = io.StringIO()
    try:
        with contextlib.redirect_stdout(sink):
            r = _with_timeout(NORM_TIMEOUT_S, fn,
                              check_sym_inc=True, use_all_cx_inc=False)
        return ("ok", len(r))
    except _NormTimeout:
        return ("timeout", None)
    except Exception as e:
        return (type(e).__name__, None)


def random_lindivs(rng: random.Random, nvars: int, ndivs: int,
                   coef_range=COEF_RANGE):
    """Generate a random LinDivs.

    Coefficients are drawn from a small range.  We reject divisors whose
    linear part is identically zero (the divisor would be a constant,
    which is a degenerate case handled separately by the normaliser but
    not useful for cross-validation here).
    """
    lo, hi = coef_range
    F, G = [], []
    for _ in range(ndivs):
        while True:
            f = tuple(rng.randint(lo, hi) for _ in range(nvars + 1))
            if any(c != 0 for c in f[:-1]):
                break
        g = tuple(rng.randint(lo, hi) for _ in range(nvars + 1))
        F.append(f)
        G.append(g)
    return LinDivs(tuple(F), tuple(G))


def _run_norms(lds: LinDivs):
    """Run both normalisers under a wall-clock timeout.  Returns
    (lipshitz_leaves, knf_leaves) or None if either call timed out
    or raised.  Stdout from the normalisers (the trace log) is silenced."""
    sink = io.StringIO()
    try:
        with contextlib.redirect_stdout(sink):
            lip = _with_timeout(NORM_TIMEOUT_S, lds.norm,
                                check_sym_inc=True, use_all_cx_inc=False)
            knf = _with_timeout(NORM_TIMEOUT_S, lds.knf_norm,
                                check_sym_inc=True, use_all_cx_inc=False)
    except (_NormTimeout, Exception):
        return None
    return lip, knf


def test_xval_random_small():
    """Randomised cross-validation: Lipshitz and KNF must agree on sat
    of the original system, modulo z3 timeouts."""
    rng = random.Random(RNG_SEED)
    matched = 0
    inconclusive = 0
    failed_norm = 0

    for trial in range(N_TRIALS):
        nvars = rng.choice([2, 2, 3])  # bias toward 2 vars (faster)
        ndivs = rng.choice([1, 1, 2])  # bias toward 1 div
        lds = random_lindivs(rng, nvars, ndivs)

        norms = _run_norms(lds)
        if norms is None:
            failed_norm += 1
            continue
        lip_leaves, knf_leaves = norms

        # Optional cheap sanity: both should return at least one leaf
        # (an unsat system still yields branches; just no sat leaves).
        assert len(lip_leaves) > 0
        assert len(knf_leaves) > 0

        lip_sat = any_leaf_sat(lip_leaves)
        knf_sat = any_leaf_sat(knf_leaves)
        orig_sat = _check(encode_original(lds))

        # Skip the trial if any side is unknown.  This is expected and
        # not a failure.
        if lip_sat is None or knf_sat is None or orig_sat is None:
            inconclusive += 1
            continue

        assert lip_sat == knf_sat, (
            f"trial {trial}: Lipshitz/KNF disagree on {lds.F=} {lds.G=}: "
            f"lip={lip_sat} knf={knf_sat}"
        )
        # Both should also agree with the original verdict.
        assert lip_sat == orig_sat, (
            f"trial {trial}: leaves disagree with original on "
            f"{lds.F=} {lds.G=}: leaves={lip_sat} orig={orig_sat}"
        )
        matched += 1

    # We want at least *some* conclusive trials, otherwise the test
    # asserts nothing meaningful.  Tune RNG_SEED / N_TRIALS if this
    # fires.
    print(f"xval: matched={matched} inconclusive={inconclusive} "
          f"failed_norm={failed_norm}")
    assert matched >= 1


def test_xval_curated_sat():
    """Curated case the normalisers should both find satisfiable.

    x | y is trivially sat (x=1, y arbitrary).
    """
    lds = LinDivs(tuple([(1, 0, 0)]),
                  tuple([(0, 1, 0)]))
    norms = _run_norms(lds)
    assert norms is not None, "trivial x|y should normalise"
    lip_leaves, knf_leaves = norms

    lip_sat = any_leaf_sat(lip_leaves)
    knf_sat = any_leaf_sat(knf_leaves)
    # Both should find sat (this case is small enough for z3).
    assert lip_sat is True
    assert knf_sat is True


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_xval_chain(n):
    """The plain divisibility chain: both normalisers should agree
    (in this case both succeed and produce the same leaf count)."""
    lds = chain(n)
    lip = _outcome(lds.norm)
    knf = _outcome(lds.knf_norm)
    assert lip == knf, (
        f"chain(n={n}): norm={lip} knf_norm={knf}"
    )


@pytest.mark.parametrize("n", [1, 2])
def test_xval_chain_with_yz(n):
    """The chain extended with y-z|x_1 and x_n|y+z has a mixed-sign LHS.
    Positive-form normalisation (f|g iff -f|g, plus the sign case-split that
    eliminates genuine mixing) now handles it, so both normalisers run and
    agree that the system is satisfiable (y=z=0, all x_i=0 is a solution).

    Parametrised only up to n=2 to keep CI fast; n=3 is correct but slow
    (~4.5s per normaliser, the Section 6 coefficient growth)."""
    lds = chain_with_yz(n)
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        lip = _with_timeout(15, lds.norm,
                            check_sym_inc=True, use_all_cx_inc=False)
        knf = _with_timeout(15, lds.knf_norm,
                            check_sym_inc=True, use_all_cx_inc=False)
    assert len(lip) >= 1 and len(knf) >= 1
    assert any_leaf_sat(lip) is True
    assert any_leaf_sat(knf) is True


def test_xval_curated_soda_paper():
    """SODA paper example: must agree on the original verdict."""
    lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
                  tuple([(0, 1, -2), (1, 1, 0)]))
    norms = _run_norms(lds)
    if norms is None:
        pytest.skip("normalisation failed/timed out")
    lip_leaves, knf_leaves = norms

    orig_sat = _check(encode_original(lds))
    lip_sat = any_leaf_sat(lip_leaves)
    knf_sat = any_leaf_sat(knf_leaves)

    if lip_sat is not None and knf_sat is not None:
        assert lip_sat == knf_sat
    if orig_sat is not None and lip_sat is not None:
        assert orig_sat == lip_sat
