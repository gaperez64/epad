import math
import itertools
from systems.linineqs import LinIneqs
from utils.matutils import Mat, Vec, vec2str, affxvars, basis_of_ker


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

    def all_ordered(self):
        for ordtyp in itertools.permutations(range(len(self.F[0]))):
            # TODO add linineqs for each subsequent pair of indices in ordtyp
            # create new lindiv and yield it, together with the order
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
            newF.append(tuple([a // d for a in f]))
            newG.append(tuple([b // d for b in g]))
        return LinDivs(tuple(newF), tuple(newG),
                       self.get_eqs(), self.get_ineqs())

    def primitive_terms(self):
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
                M = tuple([tuple(col) for col in zip(*Mt)])
                K = basis_of_ker(M)
                v[i] = math.gcd(*K[0])
            if u == tuple(v):
                return u
