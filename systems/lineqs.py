import z3

type Mat = tuple[tuple[int, ...], ...]
type Vecs = list[tuple[int, ...]]


class LinEqs:
    def __init__(self, constrs: Mat):
        self.A = constrs

    def __str__(self):
        s = ""
        for a in self.A:
            for i, c in enumerate(a):
                i += 1
                if i == len(a):
                    s += " = " + str(-1 * c)
                elif i == 1:
                    s += str(c) + "x" + str(i)
                else:
                    s += " + " + str(c) + "x" + str(i)
            s += "\n"
        return s

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
            vec = tuple([m[indets[i]] for i in range(nvars)])
            periods.append(vec)
            s.add(z3.Or([indets[i] < vec[i] for i in range(nvars)]))
        s.pop()

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
            vec = tuple([m[indets[i]] for i in range(nvars)])
            bases.append(vec)
            s.add(z3.Or([indets[i] < vec[i] for i in range(nvars)]))
        # TODO compute bases
        return bases, periods
