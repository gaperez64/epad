from systems.lineqs import LinEqs


class TestLinEqs:
    def test_str(self):
        leqs = LinEqs(tuple([(1, -1, 0, 0), (0, 1, -2, 0)]))
        assert str(leqs) ==\
            "1x1 + -1x2 + 0x3 + 0 = 0\n0x1 + 1x2 + -2x3 + 0 = 0"

    def test_solutions(self):
        leqs = LinEqs(tuple([(1, -1, 0, 0), (0, 1, -2, 0)]))
        bases, periods = leqs.solutions()
        assert len(bases) == 1
        assert len(periods) == 1
        assert bases[0] == tuple([0, 0, 0])
        assert periods[0] == tuple([2, 2, 1])
