import numpy as np
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


def affxvars(M: Mat, base: Vec, periods: Vecs) -> Mat:
    # Unpack M into numpy Ax + c
    A = np.array([row[:-1] for row in M], dtype=np.dtype(int))
    c = np.array([row[-1] for row in M], dtype=np.dtype(int))
    # Now, apply the change of variables as
    # Ny + d = A(Py + b) + c = (AP)y + (Ab + c),
    # where P is a matrix with periods as columns
    # and b is the single base
    P = np.transpose(np.array(periods, dtype=np.dtype(int)))
    b = np.array(base, dtype=np.dtype(int))
    N = A @ P
    d = A @ b + c
    # Now put back into the original data-structure types
    N = N.tolist()
    d = d.tolist()
    res = tuple([tuple(row + [d[i]]) for i, row in enumerate(N)])
    return res


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
    return B
