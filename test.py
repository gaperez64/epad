from systems.lindivs import LinDivs
from systems.lineqs import LinEqs
from systems.linineqs import LinIneqs


def main():
    lds = LinDivs(((10, 20, 30), ),
                  ((11, 22, 33), ),
                  ((1, 2, 3), (4, 5, 6)),
                  ((2, 3, 4), (5, 6, 7)))
    print(str(lds))

    leqs = LinEqs(tuple([(1, -1, 0, 0), (0, 1, -2, 0)]))
    print(str(leqs))
    print("Solutions:")
    print(leqs.solutions())

    linqs = LinIneqs(tuple([(-1, 0, -1, 10)]),
                     tuple([(1, -1, 0, 0), (0, 1, -2, 0)]))
    print(str(linqs))
    print("Solutions:")
    print(linqs.solutions())

    linqs = LinIneqs(tuple([(-1, 0, -1, 10)]))
    print(str(linqs))
    print("Solutions:")
    print(linqs.solutions())


if __name__ == "__main__":
    main()
    exit(0)
