from systems.linineqs import LinIneqs

type Mat = tuple[tuple[int, ...], ...]


class LinDivs(LinIneqs):
    def __init__(self,
                 divisors: Mat,
                 dividends: Mat,
                 eqconstrs: Mat,
                 ineqconstrs: Mat):
        assert all([len(f) == len(a) for f in divisors for a in eqconstrs])
        LinIneqs.__init__(self, eqconstrs, ineqconstrs)
        assert len(divisors) == len(dividends)
        assert all([len(f) == len(g) for f in divisors for g in dividends])
        self.F = divisors
        self.G = dividends

    def __str__(self):
        s = ""
        for (a, b) in zip(self.F, self.G):
            for i, c in enumerate(a):
                i += 1
                if i == len(a):
                    s += " + " + str(c) + " | "
                elif i == 1:
                    s += str(c) + "x" + str(i)
                else:
                    s += " + " + str(c) + "x" + str(i)
            for i, c in enumerate(b):
                i += 1
                if i == len(b):
                    s += " + " + str(c)
                elif i == 1:
                    s += str(c) + "x" + str(i)
                else:
                    s += " + " + str(c) + "x" + str(i)
            s += "\n"
        return LinIneqs.__str__(self) + s
