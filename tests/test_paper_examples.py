"""Direct integer checks for the deterministic paper examples."""

from examples.paper import exact_one_toy, fermat_ladder


def test_fermat_product_identity():
    fermat_ladder.verify(max_height=10)
    assert fermat_ladder.fermat_ratio(0) == 3
    assert fermat_ladder.fermat_ratio(1) == 5
    assert fermat_ladder.fermat_ratio(2) == 17
    assert fermat_ladder.mersenne_multiplier(3) == 255


def test_exact_one_toy_arithmetic():
    toy = exact_one_toy.build()
    assert toy.modulus == 128
    assert toy.exponent == 10
    assert toy.k == 162889462677744140625
    assert toy.c == (13, 37, -3)
    assert toy.t == (5, 109, 29)
    assert toy.rho == 107
    assert sum(row.square_condition for row in toy.rows) == 2


def test_rendering_is_deterministic():
    assert fermat_ladder.render() == fermat_ladder.render()
    assert exact_one_toy.render() == exact_one_toy.render()
