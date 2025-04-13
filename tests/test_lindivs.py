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
