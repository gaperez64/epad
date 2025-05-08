import math
from itertools import chain, combinations, permutations
from systems.linineqs import LinIneqs
from utils.matutils import (Mat, Vec,
                            vec2str, affxvars,
                            basis_of_ker, transpose)


class LinDivs(LinIneqs):
    def __init__(self,
                 divisors: Mat,
                 dividends: Mat,
                 eqconstrs: Mat = tuple(),
                 ineqconstrs: Mat = tuple()):
        assert all([len(f) == len(a) for f in divisors for a in eqconstrs])
        LinIneqs.__init__(self, eqconstrs, ineqconstrs)
        assert len(divisors) == len(dividends)
        assert all([len(f) == len(g) for f in divisors for g in dividends])
        self.F = divisors
        self.G = dividends

    def is_left_pos(self):
        return all([all([c >= 0 for c in row]) for row in self.F])

    def all_disj_left_pos(self):
        for nonneg in chain.from_iterable(combinations(self.F, n)
                                          for n in range(len(self.F) + 1)):
            ineqs = list(self.get_ineqs())
            for f in self.F:
                h = f
                if f in nonneg:
                    h = tuple([c * -1 for c in f])
                else:
                    h = list(f)
                    h[-1] += 1
                    h = tuple(h)
                assert len(h) == len(f)
                ineqs.append(h)
            ordered = LinDivs(self.F, self.G, self.get_eqs(), tuple(ineqs))
            for res in ordered.all_disj_just_divs():
                yield res.reduced()

    def all_ordered(self):
        for ordtyp in permutations(range(len(self.F[0]))):
            # TODO add linineqs for each subsequent pair of indices in ordtyp
            # create new lindiv and yield it, together with the order
            # (it could be implicit in a reordering of the variables
            # within the divisibilities)
            pass

    def get_divs(self):
        return self.F, self.G

    def __str__(self):
        divs = [vec2str(f) + " | " + vec2str(g)
                for (f, g) in zip(self.F, self.G)]
        s = LinIneqs.__str__(self)
        if len(s) != 0:
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
        print(f"order = {order}")
        not_increasing = []
        for pt in self.primitive_terms():
            assert all([c >= 0 for c in pt])
            # compute the leading variable
            lvar = -1
            larger = []
            for idx in reversed(order):
                if pt[idx] != 0:
                    lvar = idx
                else:
                    larger.append(idx)
            assert lvar >= 0
            # get the vectors in the basis of the module
            # that have no nonzero values in components larger
            # than the leading variable of the primitive term
            print(f"Primitive term = {pt}")
            print(f"Leading variable is x{lvar}")
            print(f"Larger = {larger}")
            basis = self.basis_of_divmodule(pt)
            for vec in transpose(basis):
                print(f"Is {vec} independent?")
                # Vectors with nonzero coefficients for variables larger
                # than the leading one for the primitive term can be
                # excluded
                if any([vec[i] != 0 for i in larger]):
                    print("not in intersection...")
                    continue
                # Otherwise we check for linear dependence
                fac = pt[lvar] // vec[lvar]
                if any([pt[i] * fac != vec[i]
                       for i in range(len(pt))]):
                    not_increasing.append(vec)
        return not_increasing

    def primitive_terms(self) -> Vec:
        for f in self.F:
            d = math.gcd(*f)
            yield tuple([a // d for a in f])

    def all_disj_just_divs(self):
        bases, periods = LinIneqs.solutions(self)
        for b in bases:
            F = affxvars(self.F, b, periods)
            G = affxvars(self.G, b, periods)
            yield LinDivs(F, G)

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
        # Now we construct a matrix whose colums are h and all elements of g
        # multiplied by their coefficient in v
        H = [h]
        for c, g in zip(v, self.G):
            H.append(tuple([c * a for a in g]))
        H = tuple(H)
        return transpose(H)
