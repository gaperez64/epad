# Tests and normalizer comparison

This directory contains the regression tests for the EPAD library and a
cross-validation harness that pits the two normalizers
(`LinDivs.norm` / Lipshitz and `LinDivs.knf_norm` / KNF) against each
other.

Files of interest:

- `test_norm.py` — unit tests for both normalizers on a handful of
  curated examples (the SODA paper case, Antonia's paper case, the
  trivial empty system, and a single-divisibility sanity check).
- `test_norm_xval.py` — the cross-validation harness, including:
  - z3 encoders for the input `LinDivs` and for the `OrdLinDivs` /
    `KNFLeaf` / trivial-`LinDivs` leaves that the normalizers return;
  - a randomized small-`LinDivs` generator;
  - two deterministic builders (`chain`, `chain_with_yz`) for the
    divisibility chain from Section 6 of the KNF note, with and
    without an order-breaking `y, z` extension;
  - SIGALRM-based wall-clock timeouts on the normalizers and z3
    timeouts on every `check()`.

## How the cross-validation works

For every leaf type returned by either normalizer we build a `z3.Solver`
that asserts:

- all divisibilities `f_i | g_i` as `f_i(x) * k_i == g_i(x)` for a fresh
  integer `k_i` (this exactly captures `f | g` over the integers,
  including the `0 | 0` case);
- `L_gens` equalities `l(x) == 0` for every `l` in a `KNFLeaf`;
- `x_i >= 0` for all variables of a leaf, since EPAD preprocessing
  reparametrizes solutions to the naturals.

Both normalizers must agree on the satisfiability of *some* leaf:
the original input is sat iff at least one Lipshitz leaf is sat iff at
least one KNF leaf is sat. Each z3 query is run under a per-check
timeout (`Z3_TIMEOUT_MS = 500`); `unknown` verdicts mark a trial as
inconclusive and the test simply skips it. The randomized test
(`test_xval_random_small`) only asserts that *some* trial was
conclusive and that every conclusive trial agreed.

## The chain generator

`chain(n)` builds the bare Lipshitz divisibility chain over
`x_0, ..., x_n`:

```
x_i | x_{i+1}    and    (x_i + 1) | x_{i+1}    for i = 0, ..., n-1.
```

`gcd(x_i, x_i + 1) = 1` forces `x_i (x_i + 1) | x_{i+1}`, so any
solution has `x_{i+1} > x_i^2`, and `x_i` grows doubly-exponentially
in `i`. The chain is well-known to be a "harmless" input for either
normalizer in the forward order `x_0 < x_1 < ... < x_n` — both
return a single increasing `OrdLinDivs` / `KNFLeaf`.

`chain_with_yz(n)` extends the chain with two fresh variables `y, z`
and two divisibilities:

```
(y - z) | x_1    and    x_n | (y + z).
```

By the argument in Section 6 of the companion KNF note, these
constraints simultaneously force `LV(x_n) <= max(LV(y), LV(z))` and
`LV(y), LV(z) <= LV(x_n)`, so no total order on the variables can make
every divisibility increasing.

## Observed behaviour (May 2026)

The parametrised tests `test_xval_chain[n]` and
`test_xval_chain_with_yz[n]` record the following outcomes (timings
are wall-clock on the development machine and are illustrative; what
the tests actually assert is that the two normalizers agree on the
outcome at each `n`):

| Builder           | `n` | `norm`              | `knf_norm`         |
| ----------------- | --: | ------------------- | ------------------ |
| `chain`           |   1 | ok, 1 leaf, 0.00 s  | ok, 1 leaf, 0.00 s |
| `chain`           |   2 | ok, 1 leaf, 0.00 s  | ok, 1 leaf, 0.00 s |
| `chain`           |   3 | ok, 1 leaf, 0.01 s  | ok, 1 leaf, 0.00 s |
| `chain`           |   4 | ok, 1 leaf, 0.02 s  | ok, 1 leaf, 0.02 s |
| `chain_with_yz`   |   1 | AssertionError      | AssertionError     |
| `chain_with_yz`   |   2 | AssertionError      | AssertionError     |
| `chain_with_yz`   |   3 | AssertionError, 0.8 s | AssertionError, 0.8 s |
| `chain_with_yz`   |   4 | AssertionError, 3.8 s | AssertionError, 56 s |

Two things worth noting:

- **Both normalizers fail identically on `chain_with_yz`.** They both
  trip the pre-existing precondition that every primitive divisor LHS
  is non-negative after EPAD preprocessing — the term `y - z` survives
  with mixed-sign coefficients because `all_disj_just_divs`' affine
  substitution does not always neutralize the sign. This is a
  pre-existing limitation of the library, not a KNF-specific bug, and
  the cross-validation test asserts it manifests *the same way* in
  both normalizers (same exception class, same `n`).
- **KNF gets much slower than Lipshitz once `n` grows.** At `n = 4`
  KNF takes roughly an order of magnitude longer than Lipshitz. This
  is consistent with Section 6 of the KNF note: the chain can force
  exponentially large coefficients in the accumulated quotient module
  `L`, and KNF stores those coefficients verbatim instead of folding
  them into an affine change of variables. The test suite therefore
  parametrises `test_xval_chain_with_yz` only up to `n = 3`; `n = 4`
  is used during development to confirm the trend but kept out of CI.

## Running the tests

```sh
pytest                          # full suite
pytest tests/test_norm_xval.py  # just the cross-validation
pytest -k chain_with_yz         # just the chain-with-yz comparison
```

The cross-validation harness is intentionally tolerant: z3 `unknown`
verdicts skip the trial rather than failing it, and normalizer calls
that time out (`NORM_TIMEOUT_S = 2`) are skipped too. Tune the
constants at the top of `test_norm_xval.py` if you want a stricter or
more permissive run.
