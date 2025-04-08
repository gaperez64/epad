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
