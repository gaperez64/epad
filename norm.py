from systems.lindivs import LinDivs

# x + 1 | y - 2 && x + 1 | x + y
# is NOT increasing w.r.t. x < y
order = tuple([0, 1])
lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
              tuple([(0, 1, -2), (1, 1, 0)]))
print("# INPUT SYSTEM")
print(str(lds))

print("# LEFT-POSITIVE SYSTEMS OF DIVS FORM")
for s in lds.all_disj_left_pos():
    print("## SYSTEM:")
    print(str(s))

    neqs = s.all_non_increasing(order)
    print("## ALL NON INCREASING")
    print(neqs)
