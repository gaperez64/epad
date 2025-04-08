class LinEqs:
    def __init__(self, constrs: tuple[int, ...]):
        self.A = constrs

    def __str__(self):
        s = ""
        for a in self.A:
            for i, c in enumerate(a):
                i += 1
                if c == len(a):
                    s += " = " + str(c)
                elif c == 0:
                    s += "x" + str(i)
                else:
                    s += " + x" + str(i)
        return s
