import z3
from utils.matutils import Mat, vec2str


class LinEqs:
    def __init__(self, constrs: Mat):
        self.A = constrs

    def get_eqs(self):
        return self.A

    def get_dim(self):
        assert len(self.A) > 0
        return len(self.A[0])

    def __str__(self):
        eqs = [vec2str(a) + " = 0" for a in self.A]
        return "\n".join(eqs)

    def solutions(self):
        assert len(self.A) != 0
        nvars = len(self.A[0]) - 1
        indets = [z3.Int("x" + str(i)) for i in range(nvars)]
        s = z3.Solver()
        # Make sure we focus on nonnegative values
        s.add(z3.And([indets[i] >= 0 for i in range(nvars)]))

        # First we get the periods
        s.push()
        # make sure we get nonzero solutions
        s.add(z3.Or([indets[i] >= 1 for i in range(nvars)]))
        periods = []
        needs_bases = False
        for a in self.A:
            linpoly = 0
            for i, c in enumerate(a):
                if i == len(a) - 1:
                    needs_bases = needs_bases or (c != 0)
                    s.add(linpoly == 0)
                else:
                    linpoly += str(c) * indets[i]
        while s.check() != z3.unsat:
            m = s.model()
            vec = tuple([int(m[indets[i]].as_string()) for i in range(nvars)])
            periods.append(vec)
            s.add(z3.Or([indets[i] < vec[i] for i in range(nvars)]))
        s.pop()

        # FIXME: In most applications we need the periods to be less than the
        # number of variables. This can be enforced, but for now we just check
        # that it holds
        assert len(periods) < nvars

        # Get bases
        if not needs_bases:
            return [tuple([0] * nvars)], periods
        bases = []
        for a in self.A:
            linpoly = 0
            for i, c in enumerate(a):
                if i == len(a) - 1:
                    s.add(linpoly == (-1 * c))
                else:
                    linpoly += str(c) * indets[i]
        while s.check() != z3.unsat:
            m = s.model()
            vec = tuple([int(m[indets[i]].as_string()) for i in range(nvars)])
            bases.append(vec)
            s.add(z3.Or([indets[i] < vec[i] for i in range(nvars)]))

        return bases, periods
