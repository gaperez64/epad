from systems.lindivs import LinDivs
from systems.lineqs import LinEqs


def main():
    lds = LinDivs(((1, 2, 3), (4, 5, 6)),
                  ((2, 3, 4), (5, 6, 7)),
                  ((10, 20, 30), ),
                  ((11, 22, 33), ))
    print(str(lds))

    leqs = LinEqs(tuple([(1, -1, 0, 0), (0, 1, -2, 0)]))
    print(str(leqs))
    leqs.solutions()


if __name__ == "__main__":
    main()
    exit(0)
