type Vec = tuple[int, ...]


def vec2str(v: Vec) -> str:
    s = ""
    for i, c in enumerate(v):
        i += 1
        if i == len(v):
            s += " + " + str(c)
        elif i == 1:
            s += str(c) + "x" + str(i)
        else:
            s += " + " + str(c) + "x" + str(i)
    return s
