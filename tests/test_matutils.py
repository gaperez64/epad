from systems.matutils import affxvars


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
