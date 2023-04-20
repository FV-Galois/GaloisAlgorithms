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

def deg_5(f):

    f_roots = nroots(f, n=15)

    F_G20 = lambda x: (x[0]*x[1] + x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[0] - x[0]*x[2] - x[1]*x[4] - x[4]*x[1] - x[1]*x[3] - x[3]*x[0])**2
    G20_RC = [Permutation(4), Permutation(0,1), Permutation(0,2,3,4), Permutation(0,3)(1,2,4), Permutation(1,2,3), Permutation(1,4,3,2)]
    fix_perm = compute_resolvent_roots(f_roots, F_G20, G20_RC)

    if fix_perm is None:
        if is_square(discriminant(f)):
            return 'A5'
        else:
            return 'S5'
    permute_list(f_roots, fix_perm)

    if not is_square(discriminant(f)):
        return 'G20'

    F_G5 = lambda x: x[0]*x[1]*x[1] + x[1]*x[2]*x[2] + x[2]*x[3]*x[3] + x[3]*x[4]*x[4] + x[4]*x[0]*x[0]
    G5_RC = [Permutation(4), Permutation(1,4)(2,3)]
    fix_perm = compute_resolvent_roots(f_roots, F_G5, G5_RC)

    if fix_perm is None:
        return 'AG10'
    else:
        return 'G5'

def deg_7(f):

    f_roots = nroots(f, n=15)

    F_G168 = lambda x: x[0]*x[1]*x[3] + x[0]*x[2]*x[6] + x[0]*x[4]*x[5] + x[1]*x[2]*x[4] + x[1]*x[5]*x[6] + x[2]*x[3]*x[5] + x[3]*x[4]*x[6]
    G168_RC = [Permutation(6), Permutation(0, 1), Permutation(0, 2, 3, 4, 5, 6), Permutation(0, 3, 5)(1, 2, 4, 6), Permutation(0, 2, 3, 4, 5, 6, 1), Permutation(0, 4)(1, 3, 6, 2, 5),
               Permutation(0, 3, 5, 1, 2, 4, 6), Permutation(0, 3, 5)(2, 4, 6), Permutation(0, 5, 2, 6, 3)(1, 4), Permutation(0, 2, 1), Permutation(0, 4)(1, 2, 5)(3, 6),
               Permutation(0, 3, 5, 1)(2, 4, 6), Permutation(0, 6, 4, 2)(1, 5, 3), Permutation(0, 6, 3)(1, 4)(2, 5), Permutation(0, 3, 4, 5, 6), Permutation(0, 5, 2)(1, 4)(3, 6),
               Permutation(0, 4, 1, 2, 5)(3, 6), Permutation(0, 1, 5, 3)(2, 6, 4), Permutation(1, 6, 5, 4, 3, 2), Permutation(0, 6, 4, 2, 1, 5, 3), Permutation(1, 5, 3)(2, 6, 4),
               Permutation(0, 4, 6, 1, 2, 3, 5), Permutation(0, 3, 4, 5, 6, 1), Permutation(0, 2)(1, 6, 5, 4, 3), Permutation(0, 2, 3, 4, 5, 1), Permutation(0, 1, 6, 5, 4, 3, 2),
               Permutation(0, 2, 3, 4, 5), Permutation(0, 6, 5, 1)(2, 3), Permutation(0, 3, 5, 2, 4, 6), Permutation(0, 6, 5)(2, 3)]
    fix_perm = compute_resolvent_roots(f_roots, F_G168, G168_RC)

    if fix_perm is not None:

        permute_list(f_roots, fix_perm)

        F_G21 = lambda x: x[0] * x[1] * x[3] + x[0] * x[1] * x[5] + x[0] * x[2] * x[3] + x[0] * x[2] * x[6] + x[0] * x[4] * x[5] + x[0] * x[4] * x[6] + x[1] * x[2] * x[4] + x[1] * x[2] * x[6] + x[1] * x[3] * x[4] + x[1] * x[5] * x[6] + x[2] * x[3] * x[5] + x[2] * x[4] * x[5] + x[3] * x[4] * x[6] + x[3] * x[5] * x[6]
        G21_RC = [Permutation(6), Permutation(1, 6, 3, 2)(4, 5), Permutation(1, 2, 3, 6)(4, 5), Permutation(1, 3)(2, 6), Permutation(0, 1, 3)(2, 4, 6), Permutation(1, 3)(2, 4, 6, 5), Permutation(0, 1, 4, 5, 6, 3, 2), Permutation(0, 6, 2)(1, 5, 3)]
        fix_perm = compute_resolvent_roots(f_roots, F_G21, G21_RC)
        if fix_perm is not None:

            permute_list(f_roots, fix_perm)

            F_G7 = lambda x: x[0]*x[1] + x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6] + x[6]*x[0]
            G7_RC = [Permutation(6), Permutation(1,2,4)(3,6,5), Permutation(1,4,2)(3,5,6)]
            fix_perm = (f_roots, F_G7, G7_RC)
            if fix_perm is not None:
                return 'G7'
            else:
                return 'G21'

        else:
            return 'G168'

    F_G42 = lambda x: x[0]*x[1]*x[3] + x[0]*x[1]*x[5] + x[0]*x[2]*x[3] + x[0]*x[2]*x[6] + x[0]*x[4]*x[5] + x[0]*x[4]*x[6] + x[1]*x[2]*x[4] + x[1]*x[2]*x[6] + x[1]*x[3]*x[4] + x[1]*x[5]*x[6] + x[2]*x[3]*x[5] + x[2]*x[4]*x[5] + x[3]*x[4]*x[6] + x[3]*x[5]*x[6]
    G42_RC = [Permutation(6), Permutation(0, 1), Permutation(0, 2, 3, 4, 5, 6), Permutation(0, 3, 5)(1, 2, 4, 6), Permutation(0, 2, 3, 4, 5, 6, 1),
              Permutation(0, 4)(1, 3, 6, 2, 5), Permutation(0, 3, 5, 1, 2, 4, 6), Permutation(0, 3, 5)(2, 4, 6), Permutation(0, 5, 2, 6, 3)(1, 4),
              Permutation(0, 4, 1, 3, 6, 2, 5), Permutation(0, 4)(1, 3, 6)(2, 5), Permutation(0, 4)(1, 2, 5)(3, 6), Permutation(0, 3, 5, 1)(2, 4, 6),
              Permutation(0, 6, 4, 2)(1, 5, 3), Permutation(0, 5, 2, 6, 3, 1, 4), Permutation(0, 5, 1, 4, 2, 6, 3), Permutation(0, 5, 3)(1, 4)(2, 6),
              Permutation(0, 4, 1, 3, 6)(2, 5), Permutation(0, 5, 2, 6, 4, 1, 3), Permutation(0, 4, 1, 2, 5)(3, 6), Permutation(0, 1, 5, 3)(2, 6, 4),
              Permutation(1, 6, 5, 4, 3, 2), Permutation(0, 6, 4, 2, 1, 5, 3), Permutation(0, 6, 4, 1, 5, 3, 2), Permutation(0, 6, 4, 3, 1, 5, 2),
              Permutation(0, 5)(1, 4, 2, 6, 3), Permutation(0, 4)(1, 3, 5)(2, 6), Permutation(0, 5, 3, 1, 4)(2, 6), Permutation(0, 5, 3)(1, 4, 2, 6),
              Permutation(0, 6, 5, 3, 1, 4, 2), Permutation(0, 5, 2, 6, 4)(1, 3), Permutation(0, 5, 1, 3)(2, 6, 4), Permutation(0, 2)(1, 6, 5, 4, 3),
              Permutation(1, 5, 3)(2, 6, 4), Permutation(0, 1, 6, 5, 4, 3, 2), Permutation(0, 2, 3, 4, 5), Permutation(0, 3, 5, 6, 1, 2, 4),
              Permutation(0, 6, 4)(1, 5, 3, 2), Permutation(0, 4, 6, 2, 5, 1, 3), Permutation(0, 2, 1, 6), Permutation(0, 6, 4, 3, 2)(1, 5),
              Permutation(0, 5, 2)(1, 4)(3, 6), Permutation(1, 4, 2, 6, 5, 3), Permutation(0, 3, 4, 6, 1, 5), Permutation(0, 4)(1, 2, 6, 3, 5),
              Permutation(0, 5, 3, 1, 4, 2, 6), Permutation(1, 5, 4, 3, 2), Permutation(0, 6, 5, 3)(1, 4, 2), Permutation(0, 3, 5)(1, 2, 4),
              Permutation(0, 6, 5, 2)(1, 4, 3), Permutation(0, 1, 6, 5), Permutation(0, 2, 1, 6, 5, 4, 3), Permutation(0, 2, 1), Permutation(0, 2, 3, 4, 5, 1),
              Permutation(0, 3, 5, 6)(1, 2, 4), Permutation(0, 3, 5, 6, 2, 4), Permutation(0, 4, 6, 2, 5)(1, 3), Permutation(0, 3, 4, 5, 6, 1),
              Permutation(0, 5, 3)(1, 4, 2), Permutation(1, 6, 2), Permutation(0, 5, 2, 1, 4)(3, 6), Permutation(0, 1, 5, 4, 3, 2),
              Permutation(0, 1, 5, 3, 2, 6, 4), Permutation(1, 4)(2, 6)(3, 5), Permutation(0, 6, 1, 5, 4, 3, 2), Permutation(0, 1, 6),
              Permutation(0, 2, 3, 4, 6, 5), Permutation(1, 5, 4, 2), Permutation(0, 3, 5, 1, 2, 4), Permutation(0, 4, 5, 1, 3)(2, 6),
              Permutation(0, 4)(1, 3, 6), Permutation(0, 1, 5, 6)(2, 3), Permutation(0, 3, 4, 5, 6), Permutation(0, 2), Permutation(0, 3, 5, 2, 4, 6),
              Permutation(0, 1, 5, 4, 2), Permutation(0, 5, 1, 4)(2, 6, 3), Permutation(0, 5, 3, 1, 4, 2), Permutation(0, 6, 4, 1, 5, 3),
              Permutation(1, 5, 3, 2, 6, 4), Permutation(0, 2, 4, 5)(3, 6), Permutation(0, 6)(1, 5, 4, 3, 2), Permutation(1, 6), Permutation(0, 4, 1, 2, 6, 5, 3),
              Permutation(0, 1, 6)(3, 4), Permutation(0, 3, 4, 6, 1), Permutation(0, 2, 4, 5)(1, 6), Permutation(0, 1)(2, 5, 4), Permutation(0, 4, 6, 1, 2, 3, 5),
              Permutation(0, 3, 4, 5, 6)(1, 2), Permutation(0, 2, 3, 4, 5)(1, 6), Permutation(0, 3, 5, 2, 4, 6, 1), Permutation(0, 2, 1, 6)(3, 4),
              Permutation(0, 3, 4, 6), Permutation(0, 3, 2, 1)(4, 6, 5), Permutation(1, 6, 5, 4, 2, 3), Permutation(0, 1, 6, 4, 3, 5, 2), Permutation(0, 3, 2)(5, 6),
              Permutation(0, 3, 2, 1)(5, 6), Permutation(0, 4, 1, 2, 6, 3), Permutation(1, 6, 4, 3)(2, 5), Permutation(0, 5, 1, 3, 6, 2, 4), Permutation(0, 4, 6)(1, 2, 3, 5),
              Permutation(0, 4, 6, 1, 3, 5), Permutation(0, 3, 4, 5, 6, 1, 2), Permutation(0, 2, 3, 4, 5, 1, 6), Permutation(0, 5, 4, 2)(1, 6), Permutation(0, 3, 2, 1, 5)(4, 6),
              Permutation(0, 5)(1, 3, 6, 2, 4), Permutation(0, 3, 4, 6, 5), Permutation(0, 5, 2, 4)(1, 3, 6), Permutation(0, 6, 3)(1, 4)(2, 5), Permutation(0, 5, 1, 4)(2, 3, 6),
              Permutation(0, 4, 6)(1, 3, 5), Permutation(0, 4, 5, 1, 2)(3, 6), Permutation(0, 5, 4, 2, 1, 6), Permutation(0, 5, 1, 3)(2, 4, 6), Permutation(0, 6, 2, 5, 3)(1, 4),
              Permutation(0, 6, 3, 1, 4)(2, 5), Permutation(0, 4, 6, 2)(1, 3, 5)]
    fix_perm = compute_resolvent_roots(f_roots, F_G42, G42_RC)
    if fix_perm is not None:

        permute_list(f_roots, fix_perm)

        F_G14 = lambda x: x[0]*x[1] + x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6] + x[6]*x[0]
        G14_RC = [Permutation(6), Permutation(1,3,2,6,4,5), Permutation(1,5,4,6,2,3)]
        fix_perm = compute_resolvent_roots(f_roots, F_G14, G14_RC)
        if fix_perm is not None:
            return 'G14'
        else:
            return 'G42'

    if is_square(discriminant(f)):
        return 'S7'
    else:
        return 'A7'

