type Mat = tuple[tuple[int, ...], ...]


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
