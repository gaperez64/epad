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

        if not s:  # empty system
            continue

        cxs_per_order = dict()
        inc = False
        for order in s.all_orders():
            neqs = s.all_non_increasing(order)
            if check_sym_inc and len(neqs) == 0:
                ordered.append(OrdLinDivs(s, order))
                inc = True
                break
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
