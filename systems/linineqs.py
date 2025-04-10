from systems.lineqs import LinEqs

type Mat = tuple[tuple[int, ...], ...]


class LinIneqs(LinEqs):
    def __init__(self, eqconstrs: Mat,
                 ineqconstrs: Mat):
        assert all([len(a) == len(b) for a in eqconstrs for b in ineqconstrs])
        LinEqs.__init__(self, eqconstrs)
        self.B = ineqconstrs

    def __str__(self):
        s = ""
        for a in self.B:
            for i, c in enumerate(a):
                i += 1
                if i == len(a):
                    s += " <= " + str(-1 * c)
                elif i == 1:
                    s += str(c) + "x" + str(i)
                else:
                    s += " + " + str(c) + "x" + str(i)
            s += "\n"
        return LinEqs.__str__(self) + s

    def solutions(self):
        # We first need to determine how many extra vars are needed
        xvars = len(self.B)
        assert xvars >= 1
        # Now we turn the inequalities into equalities
        eqB = []
        for i, b in enumerate(self.B):
            xtras = [0] * xvars
            xtras[i] = 1
            eqB.append(tuple(xtras + list(b)))
        # and we extend the original equalities
        exA = []
        for a in self.A:
            exA.append(tuple(([0] * xvars) + list(a)))

        # We can now make a call to get solutions of the linear system
        exsystem = LinEqs(tuple(exA + eqB))
        bases, periods = exsystem.solutions()
        # Now we have to unpack by removing the values of the extra variables
        bases = [tuple(b[xvars:]) for b in bases]
        periods = [tuple(p[xvars:]) for p in periods]

        # FIXME: In most applications we need the periods to be less than the
        # number of variables. This can be enforced, but for now we just check
        # that it holds
        assert len(periods) < len(self.B[0])

        return bases, periods
