from systems.lineqs import LinEqs
from utils.matutils import Mat, vec2str


class LinIneqs(LinEqs):
    def __init__(self,
                 ineqconstrs: Mat,
                 eqconstrs: Mat = tuple()):
        assert all([len(a) == len(b) for a in eqconstrs for b in ineqconstrs])
        LinEqs.__init__(self, eqconstrs)
        self.B = ineqconstrs

    def get_ineqs(self):
        return self.B

    def __str__(self):
        s = LinEqs.__str__(self)
        ineqs = [vec2str(b) + " <= 0" for b in self.B]
        if len(s) != 0 and len(ineqs) != 0:
            s += "\n"
        return s + "\n".join(ineqs)

    def solutions(self):
        # We first need to determine how many extra vars are needed
        xvars = len(self.B)
        # Now we turn the inequalities into equalities
        eqB = []
        for i, b in enumerate(self.B):
            xtras = [0] * xvars
            xtras[i] = 1
            eqB.append(tuple(xtras + list(b)))
        # and we extend the original equalities
        exA = []
        for a in self.get_eqs():
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
        assert len(periods) < exsystem.get_dim()

        return bases, periods
