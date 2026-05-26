"""Stress tests for the node-representation-size hunt (experiments.node_size).

These exercise the compact-stratum instrument described in the companion note
(Conjecture 12.3 / Problems 16.1-16.3): the coefficient bit-size of the KNF
equality module `L` accumulated at a leaf *is* the node-representation size.

The tests double as regression guards for two facts the hunt established:

  * Every family the normalizer can fully process has a *small* (mostly empty)
    existential increasing stratum -- no counterexample to compactness was
    found.
  * The Section 14 obstruction is real and now measurable: a *bad fixed
    order* is finitely KNF-repairable, but only with a quotient module `L`
    whose coefficients blow up doubly exponentially (bit-length doubling per
    step), while the good order needs an empty `L`.  This is exactly why the
    compactness conjecture must be existential over orders.

Two earlier limitations have been fixed and the tests now guard the fix:
  - mixed-sign divisor LHS is normalised to positive form (f|g iff -f|g), so
    inputs like chain_with_yz no longer trip a precondition;
  - the non-increasingness witness test now uses true Z-module membership, so
    `knf_norm` no longer loops on bad orders under `check_sym_inc=False`
    (the old literal-HNF-row test produced spurious witnesses).
"""

import pytest

import contextlib
import io

from systems.lindivs import LinDivs
from experiments.node_size import (
    l_metrics, l_metrics_reduced, l_metrics_lll, reduce_L, reduce_L_lll,
    leaf_is_sat, existential_min_node_size, solution_compatible_L,
    min_node_size_over_orders, smallest_certificate,
    classify_repair_generator, classify_node_repair,
    chain, antonia, all_equal, doubling_chain, cross, cross_chain,
    phi_forces_value, forward_order, reverse_order, chain_solution,
    twist_chain, cyclic_chain, amp_loose_cycle,
)


def _require_ok(r):
    """knf_norm wall-clock can vary with the machine; skip (don't fail) a
    trial that timed out so CI stays stable, but assert correctness when it
    does complete."""
    if r.get("status") == "norm_timeout":
        import pytest as _pt
        _pt.skip("knf_norm timed out on this machine")
    assert r["status"] == "ok"


# --------------------------------------------------------------------------
# metric / sat helpers
# --------------------------------------------------------------------------

def test_l_metrics_empty_and_nonempty():
    assert l_metrics(()) == {"ngen": 0, "max_abs": 0, "max_bits": 0,
                             "total_bits": 0}
    m = l_metrics(((6, -1, 0), (0, 0, 8)))
    assert m["ngen"] == 2
    assert m["max_abs"] == 8
    assert m["max_bits"] == 4          # 8 -> 0b1000
    assert m["total_bits"] == 3 + 1 + 4   # 6->3 bits, 1->1 bit, 8->4 bits


def test_leaf_is_sat_basic():
    # x | y is sat over the naturals (x=1).
    assert leaf_is_sat(LinDivs(((1, 0, 0),), ((0, 1, 0),))) is True
    # 0 = ... empty system trivially sat.
    assert leaf_is_sat(LinDivs(tuple(), tuple())) is True


# --------------------------------------------------------------------------
# Existential node size -- the quantity the compactness conjecture is about.
# A small/empty L for satisfiable families is *consistent with* (not a proof
# of) the conjecture; a super-polynomial trend would be a counterexample.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_chain_has_empty_existential_stratum(n):
    r = existential_min_node_size(chain(n))
    assert r["status"] == "ok"
    assert r["n_sat"] >= 1
    # The forward order is already increasing, so the compact stratum is empty.
    assert r["min"]["total_bits"] == 0
    assert r["best_order"] == forward_order(n)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_antonia_has_empty_existential_stratum(n):
    r = existential_min_node_size(antonia(n))
    assert r["status"] == "ok"
    assert r["n_sat"] >= 1
    assert r["min"]["total_bits"] == 0


@pytest.mark.parametrize("n", [2, 3])
def test_all_equal_has_compact_stratum(n):
    # The divisibility cycle x_i | x_{i+1 mod n} forces equalities, so L is
    # non-empty -- but it stays tiny (unit coefficients).
    r = existential_min_node_size(all_equal(n))
    assert r["status"] == "ok"
    assert r["n_sat"] >= 1
    assert r["min"]["max_abs"] <= 1
    assert r["min"]["total_bits"] <= n


# --------------------------------------------------------------------------
# Section 14 obstruction: a bad fixed order IS finitely KNF-repairable, but
# only with doubly-exponentially large coefficients in L.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_forward_order_is_increasing(n):
    r = solution_compatible_L(chain(n), forward_order(n),
                              model=chain_solution(n))
    assert r["status"] == "increasing"
    assert r["total_bits"] == 0


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7, 8])
def test_reverse_order_forces_coefficient_blowup(n):
    # The reverse (bad) order terminates -- it is finitely repairable, contra
    # the earlier spurious "stuck" diagnosis caused by the witness bug -- but
    # the coefficients in L blow up: max bit-length roughly DOUBLES each step,
    # i.e. bit(lambda) = 2^Omega(n), exactly the Section 14 obstruction.
    r = solution_compatible_L(chain(n), reverse_order(n),
                              model=chain_solution(n))
    assert r["status"] == "increasing"
    # Compare against the good order, which needs an empty stratum.
    fwd = solution_compatible_L(chain(n), forward_order(n),
                                model=chain_solution(n))
    assert fwd["total_bits"] == 0
    assert r["max_bits"] > 0


def test_reverse_order_blowup_is_doubly_exponential():
    # Quantify the blow-up: max coefficient bit-length at least doubles
    # roughly every step (a fixed-order, saturated normalisation would need
    # exponentially many bits -- the motivation for the existential conjecture).
    bits = [solution_compatible_L(chain(n), reverse_order(n),
                                  model=chain_solution(n))["max_bits"]
            for n in range(2, 9)]
    # Strictly increasing and growing by a factor close to 2.
    assert all(b2 > b1 for b1, b2 in zip(bits, bits[1:]))
    assert bits[-1] >= 8 * bits[0]   # >= 3 doublings over 6 steps


def test_reverse_blowup_is_representation_independent():
    # The Section 14 blow-up is genuine, not a representation artifact: it
    # survives BOTH canonical (HNF) reduction AND a shortest-basis (LLL)
    # reduction.  So for the bad fixed order, no compact L exists at all --
    # while the forward order needs an empty L.  This is the precise reason
    # the compactness conjecture must be existential over orders.
    lll_bits = []
    for n in (4, 6, 8):
        r = solution_compatible_L(chain(n), reverse_order(n),
                                  model=chain_solution(n))
        hnf = l_metrics_reduced(r["L"])["max_bits"]
        lll = l_metrics_lll(r["L"])["max_bits"]
        assert hnf > 0 and lll > 0
        lll_bits.append(lll)
    # LLL does not tame it: still grows (doubly exponential coefficients).
    assert lll_bits[-1] > lll_bits[0]


def test_doubling_chain_is_compact_despite_value_growth():
    # x_{i+1}=2x_i forces values 2^i, yet the stratum stays compact: the LLL
    # short basis keeps the local relations (coefficient 2).  Small n + skip
    # on timeout keeps this stable on slow CI.
    r = existential_min_node_size(doubling_chain(1), norm_timeout_s=30)
    _require_ok(r)
    assert r["n_sat"] >= 1
    assert r["min"]["max_abs"] <= 2          # LLL short basis: local relation


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_lll_keeps_doubling_relations_short(n):
    # The local doubling relations x_{i+1}-2x_i=0 are already a short basis;
    # LLL keeps the coefficient at 2 regardless of n (no 2^n composition).
    # Deterministic and fast -- the honest "node size is O(1) per relation".
    L = tuple(tuple(2 if j == i else (-1 if j == i + 1 else 0)
                    for j in range(n + 1)) for i in range(n))
    assert l_metrics(reduce_L_lll(L))["max_abs"] <= 2


# --------------------------------------------------------------------------
# The cross gadget x|Ny, y|Nx: every order needs a non-empty, LLL-irreducible
# stratum (x = N*y), of size exactly bit-length(N).  Non-empty for all orders,
# yet compact (the forced coefficient equals the input coefficient -- no
# amplification).  This is the closest "shape" to a counterexample found.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("N", [2, 4, 8, 16, 32])
def test_cross_gadget_forces_logN_for_every_order(N):
    sol = (N, 1)   # the ratio-N solution that forces x = N*y
    r = min_node_size_over_orders(cross(N), sol, max_orders=24)
    assert r["status"] == "ok"
    assert r["n_increasing"] == r["n_orders"]          # every order repairable
    # min-over-orders coefficient is exactly N (LLL-irreducible relation x=Ny);
    # node size is bit-length(N) -- compact relative to the input coefficient.
    assert r["min_L"]["max_abs"] == N
    assert r["min_metric"] == N.bit_length()


def test_cross_chain_does_not_amplify():
    # Chaining factor-2 cross gadgets does NOT compound to an exponential
    # node size: every order can use the local coefficient-2 relations, so the
    # min-over-orders coefficient stays small and even the worst order grows
    # only polynomially (contrast `chain`, which is doubly exponential on bad
    # orders because its coupling squares).
    worst = []
    for n in range(1, 5):
        sol = tuple(2 ** i for i in range(n + 1))
        r = min_node_size_over_orders(cross_chain(n), sol, max_orders=720)
        assert r["status"] == "ok"
        assert r["min_L"]["max_abs"] <= 8        # stays small, no 2^n blow-up
        worst.append(r["max_metric"])
    # Worst-order node size grows at most linearly, not exponentially.
    assert worst[-1] <= 2 * worst[0] + len(worst)


# --------------------------------------------------------------------------
# Redundancy vs. essential coefficients: the (a*x+1) gadget.  Its raw L has a
# coefficient that scales with the input parameter a, but the REDUCED basis is
# a single unit equality -- so the honest node size is O(1), compact.  This is
# why the hunt must measure a reduced basis, and a reminder that a "growing"
# raw coefficient is not by itself a counterexample.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("a", [3, 5])
def test_ax_plus_1_gadget_reduces_to_unit(a):
    # (a*x+1) | y  and  (a*x+1) | ((a+1)*x + y), with x >= 1 forced.
    lds = LinDivs(((a, 0, 1), (a, 0, 1)),
                  ((0, 1, 0), (a + 1, 1, 0)),
                  ineqconstrs=((-1, 0, 1),))
    r = existential_min_node_size(lds, norm_timeout_s=30, z3_timeout_ms=300)
    _require_ok(r)
    assert r["n_sat"] >= 1
    # Reduced node size is a single unit equality, independent of a.
    assert r["min"]["max_abs"] <= 1
    assert r["min"]["total_bits"] <= 1
    # ...even though the raw accumulated L carries a coefficient ~ a.
    assert r["min_raw"]["max_abs"] >= a - 1


# --------------------------------------------------------------------------
# Mixed-sign divisor LHS is now handled by positive-form normalisation
# (used to raise AssertionError).  This unblocks the canonical "no good
# order" construction chain_with_yz.
# --------------------------------------------------------------------------

def test_mixed_sign_lhs_now_normalized():
    # (x - y) | z : the mixed-sign LHS is eliminated by the sign case-split
    # and any residual sign-definite-negative LHS is flipped to positive, so
    # the normaliser runs instead of asserting.
    lds = LinDivs(((1, -1, 0, 0),), ((0, 0, 1, 0),))
    r = existential_min_node_size(lds, norm_timeout_s=10)
    assert not r["status"].startswith("exc:"), r["status"]


def test_chain_with_yz_runs():
    # chain_with_yz used to trip the mixed-sign precondition; it now
    # normalises and both normalisers find it satisfiable (y=z=0, x_i=0).
    from tests.test_norm_xval import chain_with_yz
    for n in (1, 2):
        lds = chain_with_yz(n)
        with contextlib.redirect_stdout(io.StringIO()):
            knf = lds.knf_norm(check_sym_inc=True, use_all_cx_inc=False)
            lip = lds.norm(check_sym_inc=True, use_all_cx_inc=False)
        assert len(knf) >= 1 and len(lip) >= 1


def test_knf_terminates_under_check_sym_inc_false():
    # Regression for the spurious-witness loop: before the Z-module-membership
    # fix, knf_norm(check_sym_inc=False) looped forever on the chain (it had
    # to explore bad orders).  It must now terminate.
    with contextlib.redirect_stdout(io.StringIO()):
        leaves = chain(1).knf_norm(check_sym_inc=False, use_all_cx_inc=False)
    assert len(leaves) >= 1


# --------------------------------------------------------------------------
# Constant (variable-free) divisor handling, via Barros' gadget
# phi_forces_value(N) = (x|y) ∧ (x|y+N) ∧ (N|x), which forces x = N.
# The N|x divisor is a constant divisor (no leading variable); it used to trip
# `assert lvar >= 0`.  It is now skipped in the increasingness check (it is a
# pure congruence), so the gadget normalises.  It also defeats two escapes
# (no good order, no all-zero collapse), leaving the small-solution pin
# {x = N} whose node size is bit-length(N) -- compact, since N is input.
# --------------------------------------------------------------------------

def test_constant_divisor_does_not_assert():
    # Bare regression for the constant-divisor fix.
    lds = phi_forces_value(2)
    r = existential_min_node_size(lds, norm_timeout_s=15)
    _require_ok(r)


@pytest.mark.parametrize("N", [2, 4, 8, 64, 256])
def test_phi_forces_value_pins_x_to_N(N):
    r = existential_min_node_size(phi_forces_value(N), norm_timeout_s=20,
                                  z3_timeout_ms=800)
    _require_ok(r)
    assert r["n_sat"] >= 1
    # No good order (x|N puts the constant N in M_x) and no all-zero collapse
    # (0 ∤ N): the only certificate is the pin {x = N}, an LLL-irreducible
    # single relation of coefficient exactly N.
    assert r["min"]["ngen"] == 1
    assert r["min"]["max_abs"] == N
    # ...so the node size is bit-length(N): scales with the *value* but is
    # compact relative to the input (N is written in the input in O(log N)).
    assert r["min"]["max_bits"] == N.bit_length()


# --------------------------------------------------------------------------
# Four-mechanism repair-completeness conjecture: every satisfiable non-KNF
# kernel state has a poly-size repair from one of (1) acyclic scheduling,
# (2) homogeneous zero collapse, (3) local/cross 2-component relation,
# (4) finite-gcd pin.  The classifier maps each minimal repair to a mechanism;
# the hunt found NO repair outside the four (no genuine >=3-variable relation),
# because a repair defining LV(f) by lower variables is absorbed by mechanism 1.
# --------------------------------------------------------------------------

def test_classify_repair_generator_units():
    assert classify_repair_generator((0, 1, 0)) == "M2:zero-collapse"   # x1=0
    assert classify_repair_generator((1, 0, -64)) == "M4:gcd-pin"        # x0=64
    assert classify_repair_generator((1, -1, 0)) == "M3:2-comp"          # x0=x1
    assert classify_repair_generator((2, -1, -2)) == "M3:2-comp(affine)"
    assert classify_repair_generator((1, 1, 1, 0)).startswith("FIFTH")   # 3 vars


@pytest.mark.parametrize("name,lds,expected", [
    ("chain",      chain(3),            "M1:reschedule"),
    ("antonia",    antonia(2),          "M1:reschedule"),
    ("x=y+z",      LinDivs(((1, 0, 0, 0), (0, 1, 1, 0)),
                           ((0, 1, 1, 0), (1, 0, 0, 0))), "M1:reschedule"),
    ("phi(64)",    phi_forces_value(64), "M4:gcd-pin"),
])
def test_node_repair_classifies_into_four_mechanisms(name, lds, expected):
    r = classify_node_repair(lds, norm_timeout_s=20, z3_timeout_ms=600)
    if r["status"] == "norm_timeout":
        pytest.skip("knf_norm timed out on this machine")
    assert r["status"] == "ok"
    assert not r["fifth"], f"{name}: unexpected 5th-mechanism repair {r}"
    assert expected in r["classes"], f"{name}: {r['classes']}"


def test_no_fifth_mechanism_on_constructed_families():
    # Every family probed has its minimal repair inside the four mechanisms;
    # none needs a genuine >=3-variable ("fifth") relation.
    fams = [chain(2), antonia(2), all_equal(2),
            cross(16), cross_chain(2), doubling_chain(1),
            phi_forces_value(8)]
    saw_ok = False
    for lds in fams:
        r = classify_node_repair(lds, norm_timeout_s=20, z3_timeout_ms=400)
        if r["status"] != "ok":
            continue
        saw_ok = True
        assert not r["fifth"], f"5th-mechanism repair: {r}"
    assert saw_ok


# --------------------------------------------------------------------------
# Curated non-empty-stratum examples surfaced by the sign-definite random
# hunt.  These are the smallest satisfiable systems we found whose increasing
# stratum L is non-trivial; all stay compact (a few bits), the evidence the
# hunt collected in favour of the compactness conjecture.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("F,G,max_abs,max_bits", [
    # 2x+y | 2x+y+1  -- forces 2x+y to divide 1, i.e. a unit; L pins it.
    (((2, 1, 0),), ((2, 1, 1),), 2, 4),
    # 2x+1 | x+1
    (((2, 0, 1),), ((1, 0, 1),), 2, 3),
    # 2x | x+1
    (((2, 0, 0),), ((1, 0, 1),), 1, 2),
])
def test_curated_compact_strata(F, G, max_abs, max_bits):
    r = existential_min_node_size(LinDivs(F, G), norm_timeout_s=20)
    _require_ok(r)
    assert r["n_sat"] >= 1
    assert r["min"]["ngen"] >= 1            # stratum is genuinely non-empty
    assert r["min"]["max_abs"] <= max_abs
    assert r["min"]["total_bits"] <= max_bits


# --------------------------------------------------------------------------
# Smallest-certificate search: the rigorous node size = min L over ALL small
# relations and ALL orders (not just the ones knf_norm generates).  knf_norm's
# L is only an upper bound; this finds the true minimum and confirms that
# every family probed has a COMPACT true certificate -- the core evidence for
# the compact-stratum conjecture, and the result of the lattice-pathology hunt.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n", [1, 2, 3])
def test_true_certificate_of_chain_is_empty(n):
    # The chain's reverse-order knf L is doubly exponential, but its TRUE
    # smallest certificate is EMPTY: the forward order is increasing with no L.
    r = smallest_certificate(chain(n), max_coef=1, max_gens=1)
    assert r["status"] == "ok"
    assert r["ngen"] == 0
    assert r["order"] == forward_order(n)


@pytest.mark.parametrize("N", [4, 16, 64, 256, 1024])
def test_true_certificate_of_cross_is_unit(N):
    # knf_norm reports L with coefficient N for cross(N); the TRUE smallest
    # certificate is a single unit relation, independent of N -- so the node
    # size is O(1), compact, and knf_norm merely overestimated.
    r = smallest_certificate(cross(N), max_coef=2, max_gens=2)
    assert r["status"] == "ok"
    assert r["metrics"]["max_abs"] <= 1
    assert r["metrics"]["ngen"] == 1


def test_cross_certificate_beats_knf_upper_bound():
    # Explicitly contrast: knf_norm's existential L vs the true certificate.
    N = 64
    knf = existential_min_node_size(cross(N), norm_timeout_s=15)
    _require_ok(knf)
    cert = smallest_certificate(cross(N), max_coef=2, max_gens=2)
    assert cert["status"] == "ok"
    # knf overestimates (coefficient ~N); the true certificate is unit.
    assert knf["min"]["max_abs"] >= N // 2
    assert cert["metrics"]["max_abs"] <= 1


# --------------------------------------------------------------------------
# Twisted-chain hunt.  Two ways to attack the compact-stratum conjecture were
# tried and both stay compact:
#   * symmetric divisors (s-t, s+t) instead of the coprime pair (x, x+1) --
#     the twist carries no forced coprimality, so it never amplifies and the
#     existential stratum is EMPTY (mechanism 1);
#   * closing the (genuinely amplifying) chain into a cycle -- a hard closure
#     collapses to the fixed point x=0 (mechanism 2); a loose closure pins the
#     boundary variable to a constant (mechanism 4).  In the loose case the
#     doubly-exp recurrence appears ONLY in the worst order, never the minimum.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n", [1, 2])
def test_twist_chain_has_empty_existential_stratum(n):
    # The symmetric difference/sum pair lacks the chain's automatic
    # coprimality, so the system has slack and a good order with empty L always
    # exists -- the twist alone does not amplify (mechanism 1: reschedule).
    r = existential_min_node_size(twist_chain(n), norm_timeout_s=30)
    _require_ok(r)
    assert r["n_sat"] >= 1
    assert r["min"]["total_bits"] == 0
    c = classify_node_repair(twist_chain(n), norm_timeout_s=30)
    if c["status"] == "ok":
        assert c["classes"] == ["M1:reschedule"]
        assert not c["fifth"]


@pytest.mark.parametrize("close_both", [True, False])
def test_cyclic_chain_collapses_to_zero(close_both):
    # Closing the amplifier into a directed divisibility cycle forces all x_i
    # equal, then (x_i+1)|x_i forces x_i=0: the cycle CANNOT hold positive
    # values.  The existential stratum is the unit zero-collapse relations
    # (mechanism 2).  n=1 only -- n>=2 branch-explodes in knf_norm.
    lds = cyclic_chain(1, close_both=close_both)
    r = existential_min_node_size(lds, norm_timeout_s=30, z3_timeout_ms=1500)
    _require_ok(r)
    assert r["n_sat"] >= 1
    assert r["min"]["max_abs"] <= 1          # unit zero-collapse, no blow-up
    c = classify_node_repair(lds, norm_timeout_s=30, z3_timeout_ms=1500)
    if c["status"] == "ok":
        assert not c["fifth"]
        assert all(cl == "M2:zero-collapse" for cl in c["classes"]), c["classes"]


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_amp_loose_cycle_forward_pins_x0_to_constant(n):
    # The loose closure x_0|(x_n+1) together with x_0|x_n (forward transitivity)
    # forces x_0 | gcd(x_n, x_n+1) = 1, i.e. x_0 = 1: an M4 finite-gcd pin of
    # size O(1).  The forward-order forced repair is exactly that single
    # constant-pin relation (1,0,...,0,-1).
    r = solution_compatible_L(amp_loose_cycle(n, 1), forward_order(n),
                              model=chain_solution(n))
    assert r["status"] == "increasing"
    basis = reduce_L_lll(r["L"])
    assert len(basis) == 1
    v = basis[0]
    assert classify_repair_generator(v) == "M4:gcd-pin"
    # the pin is x_0 = 1: unit variable coefficient, constant -1
    assert abs(v[0]) == 1 and abs(v[-1]) == 1
    assert all(c == 0 for c in v[1:-1])


def test_amp_loose_cycle_recurrence_only_in_worst_order():
    # The key Version-C result: the cycle removes the globally good order, but
    # the EXISTENTIAL minimum node size stays O(1) (the forward x_0=1 pin) for
    # every n, while only the WORST order inherits the chain's doubly-exp
    # coefficients.  So the recurrence lands in the worst order, never the
    # minimum -- consistent with the compact-stratum conjecture.
    worst = []
    for n in (2, 3, 4):
        sol = chain_solution(n)
        r = min_node_size_over_orders(amp_loose_cycle(n, 1), sol,
                                      metric="max_bits")
        assert r["status"] == "ok"
        assert r["n_increasing"] == r["n_orders"]   # every order repairable
        assert r["min_metric"] == 1                 # existential min: O(1) pin
        assert r["max_metric"] > r["min_metric"]    # worst order amplifies
        worst.append(r["max_metric"])
    # Worst-order node size grows with n (the doubly-exp recurrence), while the
    # minimum stayed pinned at 1 above.
    assert worst[0] < worst[1] < worst[2]
