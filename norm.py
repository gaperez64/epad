from systems.lindivs import LinDivs, OrdLinDivs


def norm(lds: LinDivs, check_sym_inc=True, use_all_cx_inc=True):
    print("# INPUT SYSTEM")
    print(str(lds))

    to_treat = list(lds.all_disj_left_pos())
    ordered = []
    while len(to_treat) > 0:
        s = to_treat.pop()
        print("## LEFT-POSITIVE SUBSYSTEM:")
        print(str(s))

        assert s.get_dim() > 0, "empty system?"

        cxs_per_order = dict()
        inc = False
        for order in s.all_orders():
            neqs = s.all_non_increasing(order)
            if check_sym_inc and len(neqs) == 0:
                print("### Found to be (symbolically) increasing!")
                ordered.append(OrdLinDivs(s, order))
                inc = True
                break
            cxs_per_order[order] = neqs
        if inc:
            continue

        for order, neqs in cxs_per_order.items():
            print(f"### Order: {order}")
            # If desired, we can use a single counterexample and
            # drop the rest
            if not use_all_cx_inc:
                neqs = neqs[:1]
            to_treat.extend(s.all_disj_from_noninc(neqs))
    return ordered


print("Example from SODA paper, turns out to be unsat")
# x + 1 | y - 2 && x + 1 | x + y
lds = LinDivs(tuple([(1, 0, 1), (1, 0, 1)]),
              tuple([(0, 1, -2), (1, 1, 0)]))
increasing = norm(lds, check_sym_inc=True, use_all_cx_inc=False)
for o in increasing:
    print(str(o))

print("Example from Antonia's paper")
# x, x + 1 | y && y, y + 1 | z
lds = LinDivs(tuple([(1, 0, 0, 0), (1, 0, 0, 1),
                     (0, 1, 0, 0), (0, 1, 0, 1)]),
              tuple([(0, 1, 0, 0), (0, 1, 0, 0),
                     (0, 0, 1, 0), (0, 0, 1, 0)]))
increasing = norm(lds, check_sym_inc=True, use_all_cx_inc=False)
for o in increasing:
    print(str(o))

