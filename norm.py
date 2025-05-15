from systems.lindivs import LinDivs


# x + 1 | y - 2 && x + 1 | x + y
# is NOT increasing w.r.t. x < y
order = tuple([0, 1])
lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
              tuple([(0, 1, -2), (1, 1, 0)]))
print(str(lds))
neqs = lds.all_non_increasing(order)
print(neqs)
