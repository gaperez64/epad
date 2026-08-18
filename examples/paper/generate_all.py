#!/usr/bin/env python3
"""Regenerate every deterministic LaTeX fragment used by the EPAD paper."""

from __future__ import annotations

try:
    from . import exact_one_toy, fermat_ladder
except ImportError:  # Direct execution from the repository root.
    import exact_one_toy
    import fermat_ladder


def main() -> int:
    for path in (fermat_ladder.write(), exact_one_toy.write()):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
