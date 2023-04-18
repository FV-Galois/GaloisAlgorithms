from compute_resolvent import *
from permute_list import *
from sympy import symbols, Poly, nroots, discriminant
from sympy.combinatorics import Permutation
from sympy.ntheory.primetest import is_square

def deg4(f):

    f_roots = nroots(f, n=15)

    F_G8 = lambda x: x[0]*x[2] + x[1]*x[3]
    G8_RC = [Permutation(3), Permutation(0,1), Permutation(0,2,3)]
    fix_perm = compute_resolvent_roots(f_roots, F_G8, G8_RC)
    if fix_perm is None:
        if is_square(discriminant(f)):
            return 'A4'
        else:
            return 'S4'
    permute_list(f_roots, fix_perm)

    F_G41 = lambda x: x[0]*x[1]*x[1] + x[1]*x[2]*x[2] + x[2]*x[3]*x[3] + x[3]*x[0]*x[0]
    G41_RC = [Permutation(3), Permutation(1,3)]
    fix_perm = compute_resolvent_roots(f_roots, F_G41, G41_RC)
    if fix_perm is None:
        return 'G8'
    else:
        if is_square(discriminant(f)):
            return 'AG_4^2'
        else:
            return 'G_4^1'
