import numpy as np
import math
import flint

type Vec = tuple[int, ...]
type Mat = tuple[Vec, ...]
type Vecs = list[Vec]


def vec2str(v: Vec) -> str:
    s = ""
    for i, c in enumerate(v):
        i += 1
        if i == len(v):
            s += " + " + str(c)
        elif i == 1:
            s += str(c) + "x" + str(i)
        else:
            s += " + " + str(c) + "x" + str(i)
    return s


def matmul(M: Mat, N: Mat) -> Mat:
    assert len(M[0]) == len(N)
    M = np.array(M, dtype=np.dtype(int))
    N = np.array(N, dtype=np.dtype(int))
    res = M @ N
    return tuple([tuple(row) for row in res.tolist()])


def affxvars(M: Mat, base: Vec, periods: Vecs) -> Mat:
    # Unpack M into numpy Ax + c
    A = np.array([row[:-1] for row in M], dtype=np.dtype(int))
    c = np.array([row[-1] for row in M], dtype=np.dtype(int))
    # Now, apply the change of variables as
    # Ny + d = A(Py + b) + c = (AP)y + (Ab + c),
    # where P is a matrix with periods as columns
    # and b is the single base
    b = np.array(base, dtype=np.dtype(int))
    d = A @ b + c
    if len(periods) > 0:
        P = np.transpose(np.array(periods, dtype=np.dtype(int)))
        N = A @ P
        # Now put back into the original data-structure types
        N = N.tolist()
        d = d.tolist()
        res = tuple([tuple(row + [d[i]]) for i, row in enumerate(N)])
        return res
    else:  # things are much simpler without periods
        d = d.tolist()
        res = tuple([tuple([c]) for c in d])
        return res


def transpose(M: Mat) -> Mat:
    return tuple([tuple(col) for col in zip(*M)])


def column_style_hnf(M: Mat) -> tuple[Mat, Mat]:
    A = flint.fmpz_mat(M)
    At = A.transpose()
    Ht, Tt = At.hnf(transform=True)
    H = tuple([tuple(r) for r in Ht.transpose().tolist()])
    T = tuple([tuple(r) for r in Tt.transpose().tolist()])
    return (H, T)


def basis_of_ker(M: Mat) -> Mat:
    A = flint.fmpz_mat(M)
    X, nullity = A.nullspace()
    B = tuple([tuple(r[:nullity]) for i, r in enumerate(X.tolist())])
    Bt = transpose(B)
    minB = []
    for row in Bt:
        d = math.gcd(*row)
        if d == 0:
            minB.append(row)
        else:
            minB.append(tuple([(c / math.gcd(*row)) for c in row]))
    B = transpose(tuple(minB))
    return B
