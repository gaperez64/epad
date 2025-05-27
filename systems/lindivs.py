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

    def get_dim(self):
        assert len(self.F) > 0
        return len(self.F[0]) - 1

    def is_left_pos(self):
        return all([all([c >= 0 for c in row]) for row in self.F])

    def all_disj_left_pos(self):
        fset = set(self.F)
        for nonneg in chain.from_iterable(combinations(fset, n)
                                          for n in range(len(fset) + 1)):
            ineqs = list(self.get_ineqs())
            for f in fset:
                h = f
                if f in nonneg:
                    h = tuple([c * -1 for c in f])
                else:
                    h = list(f)
                    h[-1] += 1
                    h = tuple(h)
                assert len(h) == len(f)
                ineqs.append(h)
            ordered = LinDivs(self.F, self.G, tuple(ineqs), self.get_eqs())
            for res in ordered.all_disj_just_divs():
                yield res.reduced()

    def all_disj_from_noninc(self, noninc):
        for f, g in noninc:
            S = sum([abs(c) for c in g])
            for c in range(-1 * S, S + 1):
                cf_min_g = tuple([c * f[i] - g[i]
                                  for i in range(len(f))])
                assert len(cf_min_g) == len(f)
                lds = LinDivs(self.F, self.G, self.get_ineqs(),
                              self.get_eqs() + (cf_min_g,))
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
            d = math.gcd(*f, *g)
            # we also factor out a -1 on the left if possible
            if all([c < 0 for c in f]):
                s = -1
            else:
                s = 1
            newF.append(tuple([s * a // d for a in f]))
            newG.append(tuple([b // d for b in g]))
        return LinDivs(tuple(newF), tuple(newG),
                       self.get_eqs(), self.get_ineqs())

    # The order is assumed to be increasing
    def all_non_increasing(self, order: Vec):
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
                extended.append(tuple(([0] * (idx - 1)) + [1] + ([0] * rest)))
            # Use this to extend the column-based basis of the module
            basis_of_mod = self.basis_of_divmodule(pt)
            extended = list(transpose(basis_of_mod)) + extended
            extended = transpose(tuple(extended))
            # Intersection
            ker_of_extended = basis_of_ker(extended)
            ker_of_extended = ker_of_extended[:len(basis_of_mod[0])]
            res = matmul(basis_of_mod, ker_of_extended)
            # HNF of intersection
            hnf_of_int, _ = column_style_hnf(res)
            # HNF of just pt
            hnf_of_pt, _ = column_style_hnf(tuple([tuple([c]) for c in pt]))
            hnf_of_int = transpose(hnf_of_int)
            hnf_of_pt = transpose(hnf_of_pt)
            for ridx in range(len(hnf_of_int)):
                if hnf_of_int[ridx] != hnf_of_pt[ridx]:
                    not_increasing.append(tuple([pt, hnf_of_int[ridx]]))
                    break

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
        bases, periods = LinIneqs.solutions(self)
        for b in bases:
            F = affxvars(self.F, b, periods)
            G = affxvars(self.G, b, periods)
            yield LinDivs(F, G)

    def basis_of_divmodule(self, h: Vec) -> Mat:
        v = [0] * len(self.F)
        print(f"h = {h}")
        print(str(self))
        while True:
            u = tuple(v)
            print(f"u = {u}")
            for i, f in enumerate(self.F):
                Mt = [[-1 * c for c in f], h]
                for j, g in enumerate(self.G):
                    Mt.append([u[j] * c for c in g])
                print(f"Mt = {Mt}")
                # We actually want the transpose of Mt
                M = transpose(Mt)
                print(f"M = {M}")
                K = basis_of_ker(M)
                print(f"K = {K}")
                v[i] = math.gcd(*K[0])
            if u == tuple(v):
                break
        # Now we construct a matrix whose colums are h and all elements of g
        # multiplied by their coefficient in v
        H = [h]
        for c, g in zip(v, self.G):
            H.append(tuple([c * a for a in g]))
        H = tuple(H)
        return transpose(H)
