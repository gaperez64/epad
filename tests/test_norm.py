import pytest
from systems.lindivs import LinDivs, OrdLinDivs # Ensure OrdLinDivs is imported
from norm import norm

def test_soda_paper_example():
    # Example from SODA paper, turns out to be unsat
    # x + 1 | y - 2 && x + 1 | x + y
    lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
                  tuple([(0, 1, -2), (1, 1, 0)]))
    result = norm(lds, check_sym_inc=True, use_all_cx_inc=False)
    assert result == []

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
    result = norm(lds, check_sym_inc=True, use_all_cx_inc=False)
    # This was updated based on initial test runs.
    # The original examples in norm.py didn't specify the exact output structure for this.
    # We expect a non-empty list of OrdLinDivs objects.
    assert isinstance(result, list)
    assert len(result) > 0 # Expecting one OrdLinDivs
    assert isinstance(result[0], OrdLinDivs)


def test_empty_subsystem():
    # Test case where the subsystem 's' becomes empty.
    lds = LinDivs(tuple(), tuple()) # Represents an empty system
    result = norm(lds, check_sym_inc=True, use_all_cx_inc=False)
    assert result == []

def test_soda_paper_example_all_cx():
    # Example from SODA paper, with use_all_cx_inc=True
    lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
                  tuple([(0, 1, -2), (1, 1, 0)]))
    result = norm(lds, check_sym_inc=True, use_all_cx_inc=True) # Set use_all_cx_inc to True
    assert result == []
