from systems.linineqs import LinIneqs


class TestLinIneqs:
    def test_solutions(self):
        linqs = LinIneqs(tuple([(-1, 0, -1, 10)]),
                         tuple([(1, -1, 0, 0), (0, 1, -2, 0)]))
        bases, periods = linqs.solutions()
        assert len(bases) == 1
        assert len(periods) == 1
        assert bases[0] == tuple([8, 8, 4])
        assert periods[0] == tuple([2, 2, 1])
