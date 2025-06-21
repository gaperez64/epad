from systems.lindivs import LinDivs


if __name__ == "__main__":
    # Funky example by Alessio: essentially all variables
    # should have the same value
    lds = LinDivs(tuple([(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]),
                  tuple([(0, 1, 0, 0), (0, 0, 1, 0), (1, 0, 0, 0)]))
    result = lds.norm(check_sym_inc=True, use_all_cx_inc=False)
    for r in result:
        print("--")
        print(str(r))
