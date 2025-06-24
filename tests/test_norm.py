from systems.lindivs import LinDivs, OrdLinDivs


def test_soda_paper_example():
    # Example from SODA paper, turns out to be unsat
    # x + 1 | y - 2 && x + 1 | x + y
    lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
                  tuple([(0, 1, -2), (1, 1, 0)]))
    result = lds.norm(check_sym_inc=True, use_all_cx_inc=False)
    assert len(result) > 0


def test_antonias_paper_example():
    # Example from Antonia's paper
    # x, x + 1 | y && y, y + 1 | z
    lds = LinDivs(tuple([(1, 0, 0, 0), (1, 0, 0, 1),
                         (0, 1, 0, 0), (0, 1, 0, 1)]),
                  tuple([(0, 1, 0, 0), (0, 1, 0, 0),
                         (0, 0, 1, 0), (0, 0, 1, 0)]),
                  tuple([(-1, 0, 0, 2),
                         (0, -1, 0, 2),
                         (0, 0, -1, 2)]))
    result = lds.norm(check_sym_inc=True, use_all_cx_inc=False)
    # We expect a non-empty list of OrdLinDivs objects.
    assert isinstance(result, list)
    assert len(result) > 0  # Expecting one OrdLinDivs
    assert isinstance(result[0], OrdLinDivs)


def test_empty_subsystem():
    # Test case where the subsystem 's' becomes empty.
    lds = LinDivs(tuple(), tuple())  # Represents an empty system
    result = lds.norm(check_sym_inc=True, use_all_cx_inc=False)
    assert len(result) == 1


def test_soda_paper_example_all_cx():
    # Example from SODA paper, with use_all_cx_inc=True
    lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
                  tuple([(0, 1, -2), (1, 1, 0)]))
    # Set use_all_cx_inc to True
    result = lds.norm(check_sym_inc=True, use_all_cx_inc=True)
    assert len(result) > 0
