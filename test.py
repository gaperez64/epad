from systems.lindivs import LinDivs


def main():
    lds = LinDivs(((1, 2, 3), (4, 5, 6)),
                  ((2, 3, 4), (5, 6, 7)),
                  ((10, 20, 30), ),
                  ((11, 22, 33), ))
    print(str(lds))


if __name__ == "__main__":
    main()
    exit(0)
