import math
from itertools import chain, combinations, permutations, product
from systems.linineqs import LinIneqs
from utils.matutils import (Mat, Vec, column_style_hnf,
                            vec2str, affxvars, matmul,
                            basis_of_ker, transpose)


class LinDivs(LinIneqs):
    next_id = 1

    def __init__(self,
                 divisors: Mat,
                 dividends: Mat,
                 ineqconstrs: Mat = tuple(),
                 eqconstrs: Mat = tuple()):
        assert all([len(f) == len(a) for f in divisors for a in eqconstrs])
        LinIneqs.__init__(self, ineqconstrs, eqconstrs)
        assert len(divisors) == len(dividends)
        assert all([len(f) == len(g) for f in divisors for g in dividends])
        self.F = divisors
        self.G = dividends
        self.id = LinDivs.next_id
        LinDivs.next_id += 1

    def __bool__(self):
        return len(self.F) > 0

    def get_dim(self):
        assert len(self.F) > 0
        return len(self.F[0]) - 1

    def is_left_pos(self):
        return all([all([c >= 0 for c in row]) for row in self.F])

    # This is the main function of interest in the whole codebase, it
    # implements symbolic normalization and semantic normalization; the second
    # parameter controls how aggressively we add equality constraints based on
    # witnesses of nonincreasingness
    def norm(self, check_sym_inc=True, use_all_cx_inc=True):
        log(self)
        to_treat = list(self.all_disj_left_pos())
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
                    already_ordered = OrdLinDivs(s, order)
                    ordered.append(already_ordered)
                    log(s, parent=s.id, msg=str(already_ordered))
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

    def all_disj_left_pos(self):
        some_neg = [f for f in self.F if any([c < 0 for c in f])]
        fset = set(some_neg)
        for nonneg in chain.from_iterable(combinations(fset, n)
                                          for n in range(len(fset) + 1)):
            ineqs = list(self.get_ineqs())
            for f in fset:
                h = f
                if f in nonneg:
                    h = tuple([c * -1 for c in f])
                else:
                    h = f
                ineqs.append(h)
            ordered = LinDivs(self.F, self.G, tuple(ineqs), self.get_eqs())
            log(ordered, parent=self.id, left_pos=True)
            for res in ordered.all_disj_just_divs():
                yield res.reduced()

    def all_disj_from_noninc(self, noninc):
        for f, G in noninc:
            disj_eqs = dict()
            # For each witness of nonincreasingness, we get equations that
            # encode the fact. A disjunction over subsystems per equation is
            # equisat with the original system.
            for g in G:
                disj_eqs[g] = []
                S = sum([abs(c) for c in g])
                for c in range(-1 * S, S + 1):
                    cf_min_g = tuple([c * f[i] - g[i]
                                      for i in range(len(f))])
                    assert len(cf_min_g) == len(f)
                    disj_eqs[g].append(cf_min_g)
            # For each witness, we need to choose one equation at a time to
            # get a subsytem.
            for new_eqs in product(*(disj_eqs.values())):
                lds = LinDivs(self.F, self.G, self.get_ineqs(),
                              self.get_eqs() + new_eqs)
                log(lds, parent=self.id, non_inc=True)
                for res in lds.all_disj_just_divs():
                    yield res.reduced()

    def all_orders(self):
        for ordtyp in permutations(range(len(self.F[0]) - 1)):  # skip consts
            yield ordtyp

    def get_divs(self):
        return (self.F, self.G)

    def __str__(self):
        divs = [vec2str(f) + " | " + vec2str(g)
                for (f, g) in zip(self.F, self.G)]
        s = LinIneqs.__str__(self)
        if len(s) != 0 and len(divs) != 0:
            s += "\n"
        return s + "\n".join(divs)

    def reduced(self):
        newF = []
        newG = []
        for f, g in zip(self.F, self.G):
            if all([c == 0 for c in g]):
                continue  # trivially true, everything divides 0
            d = math.gcd(*f, *g)
            assert d != 0
            # we factor out a -1 on the left if possible
            if all([c < 0 for c in f]):
                s = -1
            else:
                s = 1
            newF.append(tuple([s * a // d for a in f]))
            newG.append(tuple([b // d for b in g]))
        res = LinDivs(tuple(newF), tuple(newG),
                      self.get_eqs(), self.get_ineqs())
        log(res, parent=self.id, reduced=True)
        return res

    # The order is assumed to be increasing
    def all_non_increasing(self, order: Vec) -> Mat:
        assert len(self.F) > 0
        assert len(order) + 1 == len(self.F[0])
        not_increasing = []
        for pt in self.primitive_terms():
            assert all([c >= 0 for c in pt])
            # compute the leading variable
            lvar = -1
            lvaridx = -1
            for i, idx in enumerate(reversed(order)):
                if pt[idx] != 0:
                    lvar = idx
                    lvaridx = len(order) - 1 - i
            assert lvar >= 0
            assert 0 <= lvaridx and lvaridx < len(order)

            # NOTE: We now need a matrix representing the module spanned by
            # linear polynomials on variables smaller than or equal to the
            # leading variable. We then intersect the ones we have, and
            # compare hermite normal forms of that one and the one for
            # just pt (to check equality).

            # Matrix with row bases for all smaller variables
            extended = [tuple(([0] * len(order)) + [1])]  # and the constant(s)
            for idx in order[:lvaridx+1]:
                rest = len(order) - idx  # implicit + 1 due to constants
                assert rest + 1 + idx == len(extended[0])
                extended.append(tuple(([0] * idx) + [1] + ([0] * rest)))
            # Use this to extend the column-based basis of the module
            basis_of_mod = self.basis_of_divmodule(pt)
            extended = list(transpose(basis_of_mod)) + extended
            extended = transpose(tuple(extended))
            # Compute the intersection
            ker_of_extended = basis_of_ker(extended)
            ker_of_extended = ker_of_extended[:len(basis_of_mod[0])]
            res = matmul(basis_of_mod, ker_of_extended)
            # HNF of intersection
            hnf_of_int, a = column_style_hnf(res)
            hnf_of_int = transpose(hnf_of_int)
            # HNF of just pt
            hnf_of_pt, _ = column_style_hnf(tuple([tuple([c]) for c in pt]))
            hnf_of_pt = transpose(hnf_of_pt)

            # Get the whole set of differences
            # NOTE: the HNFs could contain trivial zeros
            non_inc_bases = [g for g in hnf_of_int
                             if g not in hnf_of_pt
                             and any([c != 0 for c in g])]
            if len(non_inc_bases) > 0:
                not_increasing.append(tuple([pt, tuple(non_inc_bases)]))

        not_increasing = tuple(not_increasing)
        return not_increasing

    def primitive_terms(self) -> Vec:
        done = set()
        for f in self.F:
            d = math.gcd(*f)
            res = tuple([a // d for a in f])
            if res not in done:
                yield res
                done.add(res)

    def all_disj_just_divs(self):
        pending = [self]
        while len(pending) > 0:
            lds = pending.pop()
            # If the system is already pure divisibilities or has no
            # divisibilities then it's ready
            if len(lds.get_eqs() +
                   lds.get_ineqs()) == 0 or not lds:
                yield lds
                continue
            # Otherwise, we need to split it
            bases, periods = LinIneqs.solutions(lds)
            for b in bases:
                F = affxvars(lds.F, b, periods)
                G = affxvars(lds.G, b, periods)
                # All LHS with constants only, become a new equality
                # constraint with an additional variable. So we count these
                # first.
                nconst = sum([1 for f in F
                              if f[-1] != 0 and all([c == 0 for c in f[:-1]])])
                eqs = []
                cleanF = []
                cleanG = []
                pref_of_zeros = tuple([0] * nconst)
                iconst = 0
                for f, g in zip(F, G):
                    if all([c == 0 for c in f]):  # f = 0
                        if all([c == 0 for c in g]):
                            # ignore it, it's trivially true
                            pass
                        else:
                            # nontrivial g and f = 0, this is only possible if
                            # g = 0 because only 0 is divisible by 0
                            eqs.append(pref_of_zeros + g)
                    elif all([c == 0 for c in f[:-1]]):  # f = c
                        pref = [0] * iconst
                        pref += [-1 * f[-1]]
                        pref += [0] * (nconst - (iconst + 1))
                        pref = tuple(pref)
                        eqs.append(pref + g)
                        # we should not be adding nonsense equations
                        assert any([c != 0 for c in eqs[-1]])
                    else:  # nonconstant f
                        cleanF.append(pref_of_zeros + f)
                        cleanG.append(pref_of_zeros + g)
                res = LinDivs(tuple(cleanF),
                              tuple(cleanG),
                              tuple(),
                              tuple(eqs))
                log(res, parent=lds.id, just_divs=True)
                pending.append(res)

    def basis_of_divmodule_knf(self, f: Vec, L_gens: Mat) -> Mat:
        """M̃f(Φ,L): least submodule of R containing {f}∪L, closed under
        divisibility transfer. Generalises basis_of_divmodule to L≠{0}."""
        v = [0] * len(self.F)
        base_gens = [f] + list(L_gens)
        while True:
            u = tuple(v)
            for i, fi in enumerate(self.F):
                # Find gcd of b s.t. b*fi ∈ span(base_gens, v[j]*g_j).
                # Encode as: -fi*b + sum(gen*λ) + sum(u[j]*g_j*z_j) = 0
                Mt = [[-1 * c for c in fi]]
                for gen in base_gens:
                    Mt.append(list(gen))
                for j, gj in enumerate(self.G):
                    Mt.append([u[j] * c for c in gj])
                M = transpose(tuple(Mt))
                K = basis_of_ker(M)
                v[i] = math.gcd(*K[0])
            if u == tuple(v):
                break
        H = list(base_gens)
        for c, g in zip(v, self.G):
            if c != 0:
                H.append(tuple([c * a for a in g]))
        return transpose(tuple(H))

    def all_non_increasing_knf(self, order: Vec, L_gens: Mat) -> Mat:
        """Return (f, witnesses) pairs where the KNF condition fails in R/L.

        KNF holds for f iff M̃f(Φ,L) ∩ (R_≤k + L) = Zf + L where k=LV_χ(f).
        """
        assert len(self.F) > 0
        assert len(order) + 1 == len(self.F[0])
        not_increasing = []
        for pt in self.primitive_terms():
            assert all([c >= 0 for c in pt])
            lvar = -1
            lvaridx = -1
            for i, idx in enumerate(reversed(order)):
                if pt[idx] != 0:
                    lvar = idx
                    lvaridx = len(order) - 1 - i
            assert lvar >= 0
            assert 0 <= lvaridx < len(order)

            # Basis of R_≤k (constant + variables appearing before lvaridx)
            r_leq_k = [tuple(([0] * len(order)) + [1])]  # constant
            for idx in order[:lvaridx + 1]:
                rest = len(order) - idx
                r_leq_k.append(tuple(([0] * idx) + [1] + ([0] * rest)))

            # M̃f(Φ,L) as column matrix
            basis_of_mod = self.basis_of_divmodule_knf(pt, L_gens)
            nbasis_cols = len(basis_of_mod[0])

            # Intersection of M̃f(Φ,L) with (R_≤k + L)
            # Combine M̃f columns and (R_≤k ∪ L) columns, find kernel
            r_leq_k_plus_L = r_leq_k + list(L_gens)
            extended = list(transpose(basis_of_mod)) + r_leq_k_plus_L
            extended = transpose(tuple(extended))
            ker = basis_of_ker(extended)
            ker_u = ker[:nbasis_cols]   # coefficients for M̃f columns
            if not any(any(c != 0 for c in row) for row in ker_u):
                continue
            res = matmul(basis_of_mod, ker_u)

            # HNF of intersection
            hnf_int, _ = column_style_hnf(res)
            hnf_int_rows = set(transpose(hnf_int))

            # HNF of Zf + L (the target)
            zf_L_cols = [pt] + list(L_gens)
            zf_L_mat = transpose(tuple(zf_L_cols))
            hnf_zf_L, _ = column_style_hnf(zf_L_mat)
            hnf_zf_L_rows = set(transpose(hnf_zf_L))

            # Witnesses: HNF basis vectors in intersection but not in Zf+L
            non_inc = [g for g in hnf_int_rows
                       if g not in hnf_zf_L_rows
                       and any(c != 0 for c in g)]
            if non_inc:
                not_increasing.append((pt, tuple(non_inc)))

        return tuple(not_increasing)

    def all_disj_from_noninc_knf(self, noninc, L_gens: Mat):
        """Yield (child_LinDivs, new_L_gens) pairs, one per case split.

        Each witness g for f generates the case split f|g, branching on
        cf - g = 0 for c ∈ {-S,...,S} with S = ‖g‖₁.  The equality is
        accumulated in L instead of being substituted into the system.
        """
        for f, G in noninc:
            disj_new_L = dict()
            for g in G:
                disj_new_L[g] = []
                S = sum(abs(c) for c in g)
                for c in range(-S, S + 1):
                    cf_min_g = tuple([c * f[i] - g[i]
                                      for i in range(len(f))])
                    disj_new_L[g].append(cf_min_g)
            for new_gens in product(*(disj_new_L.values())):
                new_L = tuple(list(L_gens) + list(new_gens))
                child = LinDivs(self.F, self.G,
                                self.get_ineqs(), self.get_eqs())
                log(child, parent=self.id, non_inc=True, L_gens=new_L)
                yield (child, new_L)

    def knf_norm(self, check_sym_inc=True, use_all_cx_inc=True):
        """KNF normalization: accumulate equalities in L instead of
        substituting variables.  Returns a list of KNFLeaf (or LinDivs
        for trivially empty branches)."""
        log(self)
        # Same sign/linear preprocessing as norm()
        to_treat = [(s, tuple()) for s in self.all_disj_left_pos()]
        ordered = []
        while to_treat:
            s, L_gens = to_treat.pop()
            if not s:
                ordered.append(s)
                continue
            cxs_per_order = dict()
            inc = False
            for order in s.all_orders():
                neqs = s.all_non_increasing_knf(order, L_gens)
                if len(neqs) == 0:
                    leaf = KNFLeaf(s, order, L_gens)
                    ordered.append(leaf)
                    log(s, parent=s.id, msg=str(leaf), L_gens=L_gens)
                    if check_sym_inc:
                        inc = True
                        break
                else:
                    cxs_per_order[order] = neqs
            if inc:
                continue
            for order, neqs in cxs_per_order.items():
                if not use_all_cx_inc:
                    neqs = neqs[:1]
                to_treat.extend(s.all_disj_from_noninc_knf(neqs, L_gens))
        return ordered

    def basis_of_divmodule(self, h: Vec) -> Mat:
        v = [0] * len(self.F)
        while True:
            u = tuple(v)
            for i, f in enumerate(self.F):
                Mt = [[-1 * c for c in f], h]
                for j, g in enumerate(self.G):
                    Mt.append([u[j] * c for c in g])
                # We actually want the transpose of Mt
                M = transpose(Mt)
                K = basis_of_ker(M)
                v[i] = math.gcd(*K[0])
            if u == tuple(v):
                break
        # Now we construct a matrix whose columns are h and all elements of g
        # multiplied by their coefficient in v
        H = [h]
        for c, g in zip(v, self.G):
            if c != 0:
                H.append(tuple([c * a for a in g]))
        H = tuple(H)
        return transpose(H)


class OrdLinDivs:
    def __init__(self,
                 divs: LinDivs,
                 order: Vec):
        assert len(divs.get_ineqs() + divs.get_eqs()) == 0
        self.F = divs.F
        self.G = divs.G
        self.order = order

    def __str__(self):
        divs = [vec2str(f) + " | " + vec2str(g)
                for (f, g) in zip(self.F, self.G)]
        o = " <= ".join([f"x{i}" for i in self.order])
        s = f"order: {o}\n"
        return s + "\n".join(divs)


class KNFLeaf:
    def __init__(self,
                 divs: LinDivs,
                 order: Vec,
                 L_gens: Mat):
        assert len(divs.get_ineqs() + divs.get_eqs()) == 0
        self.F = divs.F
        self.G = divs.G
        self.order = order
        self.L_gens = L_gens

    def __str__(self):
        divs = [vec2str(f) + " | " + vec2str(g)
                for (f, g) in zip(self.F, self.G)]
        o = " <= ".join([f"x{i}" for i in self.order])
        L_str = (", ".join(vec2str(l) + " = 0" for l in self.L_gens)
                 if self.L_gens else "empty")
        return f"order: {o}\nL: {L_str}\n" + "\n".join(divs)


log_pref = "[lindivs log] "


def log(lds: LinDivs, reduced=False, left_pos=False, non_inc=False,
        just_divs=False, parent=None, msg=None, L_gens=None):
    print(f"{log_pref}start sys {lds.id}")
    if reduced:
        print(f"{log_pref}from {parent} reduced")
    if left_pos:
        print(f"{log_pref}from {parent} left pos")
    if non_inc:
        print(f"{log_pref}from {parent} made incr.")
    if just_divs:
        print(f"{log_pref}from {parent} just divs.")
    if msg is None:
        s = str(lds)
        if L_gens:
            s += "\nL: " + ", ".join(vec2str(l) + " = 0" for l in L_gens)
        print(s)
    else:
        print(f"{log_pref}{msg}")
    print(f"{log_pref}end sys")
