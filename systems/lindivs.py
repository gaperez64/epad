import math
from itertools import chain, combinations, permutations
from systems.linineqs import LinIneqs
from utils.matutils import (Mat, Vec, column_style_hnf,
                            vec2str, affxvars, matmul,
                            basis_of_ker, transpose)


class LinDivs(LinIneqs):
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

    def __bool__(self):
        return len(self.F) > 0

    def get_dim(self):
        assert len(self.F) > 0
        return len(self.F[0]) - 1

    def is_left_pos(self):
        return all([all([c >= 0 for c in row]) for row in self.F])

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
            for res in ordered.all_disj_just_divs():
                yield res.reduced()

    def all_disj_from_noninc(self, noninc):
        for f, G in noninc:
            for g in G:
                new_eqs = []
                S = sum([abs(c) for c in g])
                for c in range(-1 * S, S + 1):
                    cf_min_g = tuple([c * f[i] - g[i]
                                      for i in range(len(f))])
                    assert len(cf_min_g) == len(f)
                    new_eqs.append(cf_min_g)

                lds = LinDivs(self.F, self.G, self.get_ineqs(),
                              self.get_eqs() + tuple(new_eqs))
                for res in lds.all_disj_just_divs():
                    yield res.reduced()

    def all_orders(self):
        for ordtyp in permutations(range(len(self.F[0]) - 1)):  # skip consts
            yield ordtyp

    def get_divs(self):
        return self.F, self.G

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
        return LinDivs(tuple(newF), tuple(newG),
                       self.get_eqs(), self.get_ineqs())

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
            hnf_of_int, _ = column_style_hnf(res)
            hnf_of_int = transpose(hnf_of_int)
            # HNF of just pt
            hnf_of_pt, _ = column_style_hnf(tuple([tuple([c]) for c in pt]))
            hnf_of_pt = transpose(hnf_of_pt)
            # Quick check: The HNFs should not contain trivial zeros
            assert all([any([c != 0 for c in g]) for g in hnf_of_pt])
            assert all([any([c != 0 for c in g]) for g in hnf_of_int])

            # Get the whole set of differences
            non_inc_bases = [g for g in hnf_of_int if g not in hnf_of_pt]
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
            if len(lds.get_eqs() + lds.get_ineqs()) == 0:
                yield lds
                continue
            # Otherwise, we need to split it
            bases, periods = LinIneqs.solutions(lds)
            for b in bases:
                F = affxvars(lds.F, b, periods)
                G = affxvars(lds.G, b, periods)
                # (1) All left-hand sides consisting of zeros only, become
                # equalities instead of divisibilities; also
                # (2) all LHS with constants only, become a new equality
                # constraint with an additional variable. So we count these
                # first.
                nconst = sum([1 for f in F if all([c == 0 for c in f[:-1]])])
                eqs = []
                cleanF = []
                cleanG = []
                pref_of_zeros = tuple([0] * nconst)
                iconst = 0
                for f, g in zip(F, G):
                    if all([c == 0 for c in f]):
                        eqs.append(pref_of_zeros + g)
                    elif all([c == 0 for c in f[:-1]]):
                        pref = [0] * iconst
                        pref += [-1 * f[-1]]
                        pref += [0] * (nconst - (iconst + 1))
                        pref = tuple(pref)
                        eqs.append(pref + g)
                    else:
                        cleanF.append(pref_of_zeros + f)
                        cleanG.append(pref_of_zeros + g)
                if len(cleanF) == 0:
                    print("We got rid of all divisibilities! solved?")
                else:
                    res = LinDivs(tuple(cleanF),
                                  tuple(cleanG),
                                  tuple(),
                                  tuple(eqs))
                    pending.append(res)

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
