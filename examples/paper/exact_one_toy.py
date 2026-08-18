#!/usr/bin/env python3
"""Generate and verify the fixed Exact-1-in-3 toy instance in the paper."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product
from math import prod
from pathlib import Path


DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[3]
    / "epad-nexp"
    / "generated"
    / "exact_one_toy.tex"
)

# A literal is (variable index, +1 for positive / -1 for negated).
CLAUSES = (
    ((0, +1), (1, +1), (2, +1)),
    ((0, -1), (1, +1), (2, -1)),
)
ELL = 3
PREFIX_LENGTH = 2
FIXED_Z = 1


@dataclass(frozen=True)
class ToyRow:
    bits: tuple[int, ...]
    signs: tuple[int, ...]
    word: int
    x0: int
    x: int
    x0_mod_m: int
    x_mod_m: int
    clauses_and_prefix: bool
    square_condition: bool


@dataclass(frozen=True)
class Toy:
    d: tuple[tuple[int, ...], ...]
    c: tuple[int, ...]
    m: int
    a: int
    modulus: int
    exponent: int
    primes: tuple[int, ...]
    k: int
    p_factors: tuple[int, ...]
    lambdas: tuple[int, ...]
    t: tuple[int, ...]
    theta: tuple[int, ...]
    h: int
    tau: int
    rho: int
    rows: tuple[ToyRow, ...]


def first_odd_primes(count: int) -> tuple[int, ...]:
    primes: list[int] = []
    candidate = 3
    while len(primes) < count:
        if all(
            candidate % divisor
            for divisor in range(3, int(candidate**0.5) + 1, 2)
        ):
            primes.append(candidate)
        candidate += 2
    return tuple(primes)


def literal_value(bits: tuple[int, ...], literal: tuple[int, int]) -> int:
    index, sign = literal
    return bits[index] if sign == 1 else 1 - bits[index]


def assignment_is_target(bits: tuple[int, ...]) -> bool:
    exact_one = all(
        sum(literal_value(bits, literal) for literal in clause) == 1
        for clause in CLAUSES
    )
    word = sum(bits[index] << index for index in range(PREFIX_LENGTH))
    return exact_one and word == FIXED_Z


def build() -> Toy:
    clause_count = len(CLAUSES)
    d = tuple(
        tuple(
            sum(sign for index, sign in clause if index == variable)
            for clause in CLAUSES
        )
        for variable in range(ELL)
    )
    a = 4**clause_count
    modulus = a * 2 ** (PREFIX_LENGTH + 1)
    c = tuple(
        sum(d[variable][clause] * 4**clause for clause in range(clause_count))
        + (a * 2**variable if variable < PREFIX_LENGTH else 0)
        for variable in range(ELL)
    )
    exponent = 2 * (clause_count + ELL)
    primes = first_odd_primes(ELL)
    k = prod(prime**exponent for prime in primes)
    p_factors = tuple(k // prime**exponent for prime in primes)
    lambdas = tuple(
        (coefficient * pow(p_factor, -1, modulus)) % modulus
        for coefficient, p_factor in zip(c, p_factors)
    )
    t = tuple(
        value if value > 0 and value % prime else value + modulus
        for value, prime in zip(lambdas, primes)
    )
    theta = tuple(value * p_factor for value, p_factor in zip(t, p_factors))
    h = sum(theta)
    geometric_sum = (a - 1) // 3
    tau = -geometric_sum + a * (
        2 * FIXED_Z - (2**PREFIX_LENGTH - 1)
    )
    rho = modulus + tau

    rows = []
    for bits in product((0, 1), repeat=ELL):
        signs = tuple(2 * bit - 1 for bit in bits)
        x0 = sum(sign * weight for sign, weight in zip(signs, theta))
        x = abs(x0)
        word = sum(bits[index] << index for index in range(PREFIX_LENGTH))
        rows.append(
            ToyRow(
                bits=bits,
                signs=signs,
                word=word,
                x0=x0,
                x=x,
                x0_mod_m=x0 % modulus,
                x_mod_m=x % modulus,
                clauses_and_prefix=assignment_is_target(bits),
                square_condition=(rho * rho - x * x) % (2 * modulus) == 0,
            )
        )

    toy = Toy(
        d=d,
        c=c,
        m=clause_count,
        a=a,
        modulus=modulus,
        exponent=exponent,
        primes=primes,
        k=k,
        p_factors=p_factors,
        lambdas=lambdas,
        t=t,
        theta=theta,
        h=h,
        tau=tau,
        rho=rho,
        rows=tuple(rows),
    )
    verify(toy)
    return toy


def verify(toy: Toy) -> None:
    for index, (coefficient, weight, prime) in enumerate(
        zip(toy.c, toy.theta, toy.primes)
    ):
        assert weight % toy.modulus == coefficient % toy.modulus
        assert weight % prime != 0
        for other_index, other_prime in enumerate(toy.primes):
            if index != other_index:
                assert weight % (other_prime**toy.exponent) == 0
    assert 2 * toy.h < toy.k

    target_signs = {
        tuple(2 * bit - 1 for bit in row.bits)
        for row in toy.rows
        if row.clauses_and_prefix
    }
    assert target_signs
    for row in toy.rows:
        assert row.x <= toy.h
        assert (toy.h * toy.h - row.x * row.x) % toy.k == 0
        # Taking |X_0| identifies alpha and -alpha.  The square condition is
        # therefore true exactly when one orientation is a target assignment.
        assert row.square_condition == (
            row.signs in target_signs
            or tuple(-sign for sign in row.signs) in target_signs
        )


def sign_string(signs: tuple[int, ...]) -> str:
    return "(" + ",".join("+" if sign > 0 else "-" for sign in signs) + ")"


def render() -> str:
    toy = build()
    table_rows = []
    for row in toy.rows:
        bits = "".join(str(bit) for bit in row.bits)
        mark = "\\checkmark" if row.square_condition else "\\times"
        table_rows.append(
            "  $"
            + bits
            + "$ & $"
            + sign_string(row.signs)
            + "$ & $"
            + str(row.word)
            + "$ & $"
            + str(row.x0_mod_m)
            + "$ & $"
            + str(row.x_mod_m)
            + "$ & $"
            + mark
            + "$ \\\\"
        )

    return "\n".join(
        [
            "% Generated by ../epad-normalizer/examples/paper/exact_one_toy.py.",
            "% Do not edit by hand; regenerate with generate_all.py.",
            "For this instance,",
            "\\[",
            f"  M={toy.modulus},\\quad e={toy.exponent},\\quad",
            f"  K=(3\\cdot5\\cdot7)^{{10}}={toy.k},\\quad",
            f"  \\rho(1)={toy.rho},",
            "\\]",
            "and, writing \\(P_i=K/p_i^{10}\\),",
            "\\[",
            f"  (\\theta_0,\\theta_1,\\theta_2)=({toy.t[0]}P_0,{toy.t[1]}P_1,{toy.t[2]}P_2),",
            f"  \\qquad H={toy.h}.",
            "\\]",
            "\\begin{table}[H]",
            "\\centering",
            "\\small",
            "\\setlength{\\tabcolsep}{3.5pt}",
            "\\begin{tabular}{cccccc}",
            "\\toprule",
            "\\(u_0u_1u_2\\) & \\(\\vec\\alpha\\) & \\(\\word_2(\\vec u)\\) &",
            "\\(X_0\\bmod M\\) & \\(X\\bmod M\\) & \\(2M\\mid\\rho(1)^2-X^2\\) \\\\",
            "\\midrule",
            *table_rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\caption{All sign vectors for the fixed toy instance",
            "\\(C_0=(u_0,u_1,u_2)\\),",
            "\\(C_1=(\\neg u_0,u_1,\\neg u_2)\\), and \\(Z=1\\).",
            "Here \\(X_0=\\sum_i\\alpha_i\\theta_i\\) and \\(X=|X_0|\\).",
            "The two checked rows are complementary and yield the same \\(X\\);",
            "the proof chooses the orientation congruent to \\(\\rho(1)\\).}",
            "\\label{tab:exact-one-toy}",
            "\\end{table}",
            "",
        ]
    )


def write(output: Path = DEFAULT_OUTPUT) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render(), encoding="utf-8")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(write(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
