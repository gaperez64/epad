from systems.linineqs import LinIneqs
from systems.lindivs import LinDivs


class TestLinDivs:
    def test_str(self):
        lds = LinDivs(((10, 20, 30), ),
                      ((11, 22, 33), ),
                      ((1, 2, 3), (4, 5, 6)),
                      ((2, 3, 4), (5, 6, 7)))
        assert str(lds) ==\
            "2x1 + 3x2 + 4 = 0\n5x1 + 6x2 + 7 = 0\n"\
            "1x1 + 2x2 + 3 <= 0\n4x1 + 5x2 + 6 <= 0\n"\
            "10x1 + 20x2 + 30 | 11x1 + 22x2 + 33"

    def test_reduced(self):
        lds = LinDivs(tuple([(20, 40, 60)]),
                      tuple([(22, 44, 66)]))
        red = lds.reduced()
        f, g = red.get_divs()
        assert len(f) == 1
        assert len(g) == 1
        assert f[0] == (10, 20, 30)
        assert g[0] == (11, 22, 33)

        # Check that the -1 is also factored out
        lds = LinDivs(tuple([(-20, -40, -60)]),
                      tuple([(22, 44, 66)]))
        red = lds.reduced()
        f, g = red.get_divs()
        assert len(f) == 1
        assert len(g) == 1
        assert f[0] == (10, 20, 30)
        assert g[0] == (11, 22, 33)

        # A more elaborate test
        lds = LinDivs(tuple([(20, 40, 60),
                             (7, 14, 21)]),
                      tuple([(22, 44, 66),
                             (28, 35, 42)]))
        red = lds.reduced()
        f, g = red.get_divs()
        assert len(f) == 2
        assert len(g) == 2
        assert f[0] == (10, 20, 30)
        assert g[0] == (11, 22, 33)
        assert f[1] == (1, 2, 3)
        assert g[1] == (4, 5, 6)

    def test_primitive_terms(self):
        lds = LinDivs(tuple([(20, 40, 60), (33, 27, 99)]),
                      tuple([(22, 44, 66), (5, 7, 11)]))
        pts = [t for t in lds.primitive_terms()]
        assert len(pts) == 2
        assert (1, 2, 3) in pts
        assert (11, 9, 33) in pts

    def test_all_disj_just_divs(self):
        linqs = LinIneqs(tuple([(-1, -1, -1, 10)]),
                         tuple([(1, 1, 0, -1)]))
        lds = LinDivs(tuple([(5, 7, 0, 0)]),
                      tuple([(0, 11, -3, 1)]),
                      linqs.get_ineqs(),
                      linqs.get_eqs())
        disj = [d.get_divs() for d in lds.all_disj_just_divs()]
        assert all([len(f) == 1 for (f, g) in disj])
        disj = [(f[0], g[0]) for (f, g) in disj]
        assert ((0, 5), (-3, -26)) in disj
        assert ((0, 7), (-3, -15)) in disj

    def test_all_disj_left_pos(self):
        linqs = LinIneqs(tuple([(-1, -1, -1, 10)]),
                         tuple([(1, 1, 0, -1)]))
        lds = LinDivs(tuple([(5, 7, 0, 0)]),
                      tuple([(0, 11, -3, 1)]),
                      linqs.get_ineqs(),
                      linqs.get_eqs())
        disj = [d for d in lds.all_disj_left_pos()]
        assert all([d.is_left_pos() for d in disj])

    def test_all_non_increasing(self):
        # Examples from the introduction of our SODA paper
        # x + 1 | y - 2 is increasing w.r.t. x < y
        lds = LinDivs(tuple([(1, 0, 1)]),
                      tuple([(0, 1, -2)]))
        order = tuple([0, 1])
        assert 0 == len(lds.all_non_increasing(order))

        # x + 1 | y - 2 && x + 1 | x + y
        # is NOT increasing w.r.t. x < y
        lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
                      tuple([(0, 1, -2), (1, 1, 0)]))
        assert 0 < len(lds.all_non_increasing(order))
