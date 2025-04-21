from utils.matutils import affxvars, column_style_hnf
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
