from utils.matutils import (affxvars, column_style_hnf, basis_of_ker,
                            in_module, matmul)
import flint


class TestMatutils:
    def test_affxvars(self):
        bases = [tuple([0, 0, 0])]
        periods = [tuple([2, 2, 1])]
        r = affxvars(tuple([tuple([1, 1, 0, 0])]),
                     bases[0], periods)
        assert len(r) == 1
        assert r[0] == tuple([4, 0])
        r = affxvars(tuple([tuple([1, 0, 1, 0])]),
                     bases[0], periods)
        assert len(r) == 1
        assert r[0] == tuple([3, 0])

    def test_column_style_hnf(self):
        Al = [tuple([4, 0, 0]),
              tuple([0, 1, 0]),
              tuple([-4, 2, 0]),
              tuple([-8, 3, 0])]
        A = tuple(Al)
        H, T = column_style_hnf(A)
        assert H == A
        assert flint.fmpz_mat(H) == flint.fmpz_mat(A) * flint.fmpz_mat(T)

    def test_basis_of_ker(self):
        Al = [tuple([1, 2, 3, 4, 5]),
              tuple([6, 7, 8, 9, 10]),
              tuple([11, 12, 13, 14, 15])]
        A = tuple(Al)
        B = basis_of_ker(A)
        assert len(B) == 5
        assert len(B[0]) == 3
        exl = [tuple([1, 2, 3]),
               tuple([-2, -3, -4]),
               tuple([1, 0, 0]),
               tuple([0, 1, 0]),
               tuple([0, 0, 1])]
        ex = tuple(exl)
        assert B == ex

    def test_basis_of_ker_empty_kernel_preserves_shape(self):
        # 3x3 full-rank matrix -> trivial kernel (nullity 0).  The result
        # must still report cols(M) "rows" (one empty tuple per column) so
        # that downstream `matmul(M', basis_of_ker(M))` can rely on
        # len(basis_of_ker(M)) == cols(M).  Returning a bare () here would
        # break LinDivs.all_non_increasing, which slices and matmuls the
        # kernel against basis_of_divmodule.
        A = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        B = basis_of_ker(A)
        assert len(B) == 3, "must have one row per column of M"
        for row in B:
            assert row == (), "each row is the empty zero-column basis"
        # Downstream pipeline must accept this empty result without errors.
        Mt = ((1, 2, 3), (4, 5, 6))   # shape (2, 3)
        prod = matmul(Mt, B)
        assert len(prod) == 2 and all(r == () for r in prod)
        # column_style_hnf must also accept (m, 0) input.
        H, _T = column_style_hnf(prod)
        assert all(r == () for r in H)

    def test_basis_of_ker_zero_columns_returns_empty(self):
        # M with shape (k, 0): no columns, so the kernel space is Z^0 = {0}
        # and the result is the empty tuple.
        M = ((), (), ())
        B = basis_of_ker(M)
        assert B == ()

    def test_in_module(self):
        # Module spanned by (6,-1,0): exactly the integer multiples.
        gens = ((6, -1, 0),)
        assert in_module(gens, (6, -1, 0)) is True
        assert in_module(gens, (12, -2, 0)) is True      # 2*(6,-1,0)
        assert in_module(gens, (0, 0, 0)) is True        # zero always in
        assert in_module(gens, (6, 1, 0)) is False
        assert in_module(gens, (1, 0, 0)) is False
        # The spurious-witness case that made KNF loop: g=(6,1,12) is in
        # Zf+L for f=(1,0,1), L=⟨(6,-1,0)⟩, even though it is not an HNF row.
        zf_L = ((1, 0, 1), (6, -1, 0))
        assert in_module(zf_L, (6, 1, 12)) is True
        assert in_module(zf_L, (0, 2, 12)) is True
        # Empty module contains only the zero vector.
        assert in_module((), (0, 0)) is True
        assert in_module((), (1, 0)) is False

    def test_matmul_bignum(self):
        # Object dtype: must not overflow int64 (KNF accumulates huge coeffs).
        big = 10 ** 40
        r = matmul(((big, 0), (0, big)), ((1,), (1,)))
        assert r == ((big,), (big,))
