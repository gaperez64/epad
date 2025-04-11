from systems.linineqs import LinIneqs
from systems.matutils import vec2str

type Mat = tuple[tuple[int, ...], ...]


class LinDivs(LinIneqs):
    def __init__(self,
                 divisors: Mat,
                 dividends: Mat,
                 eqconstrs: Mat = tuple(),
                 ineqconstrs: Mat = tuple()):
        assert all([len(f) == len(a) for f in divisors for a in eqconstrs])
        LinIneqs.__init__(self, eqconstrs, ineqconstrs)
        assert len(divisors) == len(dividends)
        assert all([len(f) == len(g) for f in divisors for g in dividends])
        self.F = divisors
        self.G = dividends

    def __str__(self):
        divs = [vec2str(f) + " | " + vec2str(g)
                for (f, g) in zip(self.F, self.G)]
        s = LinIneqs.__str__(self)
        if len(s) != 0:
            s += "\n"
        return s + "\n".join(divs)
