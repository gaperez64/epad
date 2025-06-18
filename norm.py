from systems.lindivs import LinDivs, OrdLinDivs


# This is the main function of interest in the whole codebase, it implements
# symbolic normalization and semantic normalization; the second parameter
# controls how aggressively we add equality constraints based on witnesses of
# nonincreasingness
def norm(lds: LinDivs, check_sym_inc=True, use_all_cx_inc=True):
    to_treat = list(lds.all_disj_left_pos())
    ordered = []
    while len(to_treat) > 0:
        s = to_treat.pop()

        if not s:  # trivial system: no divs or just divs
            ordered.append(s)
            continue

        cxs_per_order = dict()
        inc = False
        for order in s.all_orders():
            neqs = s.all_non_increasing(order)
            if len(neqs) == 0:
                ordered.append(OrdLinDivs(s, order))
                # If, additionally, we only care about symbolic
                # increasingness, we can stop here!
                if check_sym_inc:
                    inc = True
                    break
            else:
                cxs_per_order[order] = neqs
        if inc:
            continue

        for order, neqs in cxs_per_order.items():
            # If desired, we can use a single counterexample and
            # drop the rest
            if not use_all_cx_inc:
                neqs = neqs[:1]
            to_treat.extend(s.all_disj_from_noninc(neqs))
    return ordered


if __name__ == "__main__":
    # Funky example by Alessio: essentially all variables
    # should have the same value
    lds = LinDivs(tuple([(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]),
                  tuple([(0, 1, 0, 0), (0, 0, 1, 0), (1, 0, 0, 0)]))
    print(str(lds))
    result = norm(lds, check_sym_inc=True, use_all_cx_inc=False)
    for r in result:
        print("--")
        print(str(r))
