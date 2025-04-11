from systems.matutils import affxvars
from systems.lineqs import LinEqs
from systems.lindivs import LinDivs


class TestMatutils:
    def test_vec2str(self):
        leqs = LinEqs(tuple([(1, -1, 0, 0), (0, 1, -2, 0)]))
        assert str(leqs) ==\
            "1x1 + -1x2 + 0x3 + 0 = 0\n0x1 + 1x2 + -2x3 + 0 = 0"
        lds = LinDivs(((10, 20, 30), ),
                      ((11, 22, 33), ),
                      ((1, 2, 3), (4, 5, 6)),
                      ((2, 3, 4), (5, 6, 7)))
        assert str(lds) ==\
            "2x1 + 3x2 + 4 = 0\n5x1 + 6x2 + 7 = 0\n"\
            "1x1 + 2x2 + 3 <= 0\n4x1 + 5x2 + 6 <= 0\n"\
            "10x1 + 20x2 + 30 | 11x1 + 22x2 + 33"

    def test_affxvars(self):
        pass
