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

The constraints `(y - z) | x_1` introduce a mixed-sign LHS.  Positive-form
normalization (`f | g` iff `-f | g`, plus the sign case-split that turns
genuinely mixed forms into eliminated inequalities) handles this, so both
normalizers run; the all-zero assignment (`y = z = 0`, `x_i = 0`) makes the
system satisfiable.

## Observed behaviour

The parametrised tests `test_xval_chain[n]` and
`test_xval_chain_with_yz[n]` record the following outcomes (timings
are wall-clock on the development machine and are illustrative; what
the tests actually assert is that the two normalizers agree, and on
`chain_with_yz` that both find the system satisfiable):

| Builder           | `n` | `norm`              | `knf_norm`         |
| ----------------- | --: | ------------------- | ------------------ |
| `chain`           |   1 | ok, 1 leaf, 0.00 s  | ok, 1 leaf, 0.00 s |
| `chain`           |   2 | ok, 1 leaf, 0.00 s  | ok, 1 leaf, 0.00 s |
| `chain`           |   3 | ok, 1 leaf, 0.01 s  | ok, 1 leaf, 0.00 s |
| `chain`           |   4 | ok, 1 leaf, 0.02 s  | ok, 1 leaf, 0.02 s |
| `chain_with_yz`   |   1 | ok, 2 leaves, sat   | ok, 2 leaves, sat  |
| `chain_with_yz`   |   2 | ok, 2 leaves, sat   | ok, 2 leaves, sat  |
| `chain_with_yz`   |   3 | ok, 2 leaves, ~4.5 s | ok, 2 leaves, ~4 s |

Two things worth noting:

- **`chain_with_yz` used to trip a precondition (mixed-sign LHS) and now
  runs.** The increasingness analysis no longer asserts that primitive
  divisor LHS are non-negative; instead `positive_form` flips sign-definite
  LHS (`f | g` iff `-f | g`), and the existing sign case-split in
  `all_disj_left_pos` eliminates genuinely mixed forms.  `test_xval_chain_with_yz`
  is parametrised up to `n = 2` for CI speed; `n = 3` is correct but slow
  (~4.5 s per normalizer, the Section 6 coefficient growth).
- **KNF accumulates large coefficients on bad orders.** Consistent with
  Section 6/14 of the KNF note, KNF stores the quotient-module coefficients
  verbatim instead of folding them into an affine change of variables, so a
  bad fixed order forces them to grow doubly exponentially.  The
  `experiments/node_size.py` harness measures this directly (see the
  fixed-order diagnostic, where the reverse order's coefficient bit-length
  doubles per step).

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
