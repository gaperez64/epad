from systems.lindivs import LinDivs, OrdLinDivs

# x + 1 | y - 2 && x + 1 | x + y
lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
              tuple([(0, 1, -2), (1, 1, 0)]))
print("# INPUT SYSTEM")
print(str(lds))

print("# LEFT-POSITIVE SYSTEMS OF DIVS FORM")
to_treat = list(lds.all_disj_left_pos())
ordered = []
while len(to_treat) > 0:
    s = to_treat.pop()
    print("## SYSTEM:")
    print(str(s))

    if s.get_dim() == 0:
        ordered.append(tuple([s, ()]))
        continue

    for order in s.all_orders():
        print(f"### Order: {order}")
        neqs = s.all_non_increasing(order)
        if len(neqs) > 0:
            to_treat.extend(s.all_disj_from_noninc(neqs))
        else:
            print("### Increasing!")
            ordered.append(OrdLinDivs(s, order))

for o in ordered:
    print(str(o))
