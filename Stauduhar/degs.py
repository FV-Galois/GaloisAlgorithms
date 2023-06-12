from compute_resolvent import *
from permute_list import *
from sympy import symbols, Poly, nroots, discriminant
from sympy.combinatorics import Permutation
from sympy.ntheory.primetest import is_square

def deg_3(f):
    if is_square(discriminant(f)):
        return 'A3'
    else:
        return 'S3'
def deg_4(f):

    f_roots = f.all_roots()

    F_G8 = lambda x: x[0]*x[2] + x[1]*x[3]
    G8_RC = [Permutation(3), Permutation(0,1), Permutation(2,3)]
    fix_perm = return_integer_permutation(f_roots, F_G8, G8_RC)
    if fix_perm is None:
        if is_square(discriminant(f)):
            return 'A4'
        else:
            return 'S4'
    f_roots = permute_list(f_roots, fix_perm)

    F_G41 = lambda x: x[0]*x[1]*x[1] + x[1]*x[2]*x[2] + x[2]*x[3]*x[3] + x[3]*x[0]*x[0]
    G41_RC = [Permutation(3), Permutation(0,1)(2,3)]
    fix_perm = return_integer_permutation(f_roots, F_G41, G41_RC)
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
    G20_RC = [Permutation(4), Permutation(0,1)(2,3), Permutation(0,1,3,2,4), Permutation(0,4,1,3,2), Permutation(0,1,3,4,2), Permutation(0,1,4,3,2)]
    fix_perm = return_integer_permutation(f_roots, F_G20, G20_RC)

    if fix_perm is None:
        if is_square(discriminant(f)):
            return 'A5'
        else:
            return 'S5'
    f_roots = permute_list(f_roots, fix_perm)

    if not is_square(discriminant(f)):
        return 'G20'

    F_G5 = lambda x: x[0]*x[1]*x[1] + x[1]*x[2]*x[2] + x[2]*x[3]*x[3] + x[3]*x[4]*x[4] + x[4]*x[0]*x[0]
    G5_RC = [Permutation(4), Permutation(0,1)(2,4)]
    fix_perm = return_integer_permutation(f_roots, F_G5, G5_RC)

    if fix_perm is None:
        return 'AG10'
    else:
        return 'G5'

def deg_6(f):

    f_roots = nroots(f, n=15)

    F_G120 = lambda x: (x[0]*x[1]+x[2]*x[4]+x[3]*x[5])*(x[0]*x[2]+x[3]*x[4]+x[1]*x[5])*(x[2]*x[3]+x[0]*x[5]+x[1]*x[4])*(x[0]*x[4]+x[1]*x[3]+x[2]*x[5])*(x[0]*x[3]+x[1]*x[2]+x[4]*x[5])
    G120_RC = [Permutation(5), Permutation(0,2), Permutation(1,2), Permutation(0,1,2), Permutation(0,2,1), Permutation(0,1)]
    fix_perm = return_integer_permutation(f_roots, F_G120, G120_RC)

    if fix_perm is not None:
        if is_square(discriminant(f)):
            return 'AG60'
        return 'G120'

    F_G72 = lambda x: x[0]*x[1]*x[2]+x[3]*x[4]*x[5]
    G72_RC = [Permutation(5), Permutation(1,4,3,2), Permutation(1,2,5)(3,4), Permutation(1,4,3,2,5), Permutation(1,4)(2,3), Permutation(1,3,4,2), Permutation(1,4),
              Permutation(1,2,3,4), Permutation(1,3,4,2,5), Permutation(2,5,3,4)]
    fix_perm = return_integer_permutation(f_roots, F_G72, G72_RC)

    if fix_perm is not None:

        f_roots = permute_list(f_roots, fix_perm)

        F_G236 = lambda x: (x[0]-x[1])*(x[1]-x[2])*(x[2]-x[0])*(x[3]-x[4])*(x[4]-x[5])*(x[5]-x[3])
        G236_RC = [Permutation(5), Permutation(4,5)]
        fix_perm = return_integer_permutation(f_roots, F_G236, G236_RC)

        if fix_perm is not None:

            f_roots = permute_list(f_roots, fix_perm)

            F_G18 = lambda x: (x[0]-x[1])*(x[1]-x[2])*(x[2]-x[0])+(x[3]-x[4])*(x[4]-x[5])*(x[5]-x[3])
            G18_RC = [Permutation(5), Permutation(0,1)(3,4), Permutation(4,5), Permutation(0,1)(3,5,4)]
            fix_perm = return_integer_permutation(f_roots, F_G18, G18_RC)

            if fix_perm is not None:

                f_roots = permute_list(f_roots, fix_perm)

                F_G16 = lambda x: x[0]*x[3]+x[1]*x[5]+x[2]*x[4]
                G16_RC = [Permutation(5), Permutation(0,1,2), Permutation(0,2,1)]
                fix_perm = return_integer_permutation(f_roots, F_G16, G16_RC)

                if fix_perm is not None:
                    return 'G^1_6'

                F_G26 = lambda x: x[0]*x[5]*x[5]+x[1]*x[3]*x[3]+x[2]*x[4]*x[4]+x[3]*x[1]*x[1]+x[4]*x[0]*x[0]+x[5]*x[1]*x[1]
                G26_RC = [Permutation(5), Permutation(0, 1, 2), Permutation(0, 2, 1)]
                fix_perm = return_integer_permutation(f_roots, F_G26, G26_RC)

                if fix_perm is not None:
                    return 'G^2_6'
                else:
                    return 'G18'

            else:
                if is_square(discriminant(f)):
                    return 'G^2_36'
                return 'G^1_12'
        else:
            if is_square(discriminant(f)):
                return 'AG^1_36'
            return 'G72'

    F_G48 = lambda x: x[0]*x[1]+x[2]*x[3]+x[4]*x[5]
    G48_RC = [Permutation(5), Permutation(1,3,5,2,4), Permutation(1,5)(2,4), Permutation(2,4,3), Permutation(1,2,3,4), Permutation(1,4,2), Permutation(2,3,4), Permutation(1,4,5)(2,3),
              Permutation(1,5,3,2,4), Permutation(1,2,3,5), Permutation(1,2,3), Permutation(1,4)(2,5), Permutation(1,3,2,4), Permutation(1,3)(2,4), Permutation(1,5,4,3,2)]
    fix_perm = return_integer_permutation(f_roots, F_G48, G48_RC)

    if fix_perm is not None:

        f_roots = permute_list(f_roots, fix_perm)

        F_G124 = lambda x: (x[0]+x[1]-x[2]-x[3])*(x[1]+x[3]-x[4]-x[5])*(x[4]-x[5]-x[0]-x[5])*(x[0]-x[1])*(x[2]-x[3])*(x[4]-x[5])
        G124_RC = [Permutation(5), Permutation(0,1)]
        fix_perm = return_integer_permutation(f_roots, F_G124, G124_RC)

        if fix_perm is not None:
            return 'G^1_24'

        F_G224 = lambda x: (x[0]+x[1]-x[2]-x[3])*(x[1]+x[3]-x[4]-x[5])*(x[4]+x[5]-x[0]-x[1])
        G224_RC = [Permutation(5), Permutation(0,2)(1,3)]
        fix_perm = return_integer_permutation(f_roots, F_G224, G224_RC)

        if fix_perm is not None:
            return 'G^2_24'

        if is_square(discriminant(f)):

            F_G212 = lambda x:  (x[0]+x[1]-x[2]-x[3])*(x[1]+x[3]-x[4]-x[5])*(x[4]+x[5]-x[0]-x[1])
            G212_RC = [Permutation(5), Permutation(0,2)(1,3)]
            fix_perm = return_integer_permutation(f_roots, F_G212, G212_RC)

            if fix_perm is not None:
                return 'G^2_12'
            else:
                return 'AG^3_24'
        else:
            return 'G48'

    if is_square(discriminant(f)):
        return 'A6'
    return 'S6'

def deg_7(f):

    f_roots = nroots(f, n=15)

    F_G168 = lambda x: x[0]*x[1]*x[3] + x[0]*x[2]*x[6] + x[0]*x[4]*x[5] + x[1]*x[2]*x[4] + x[1]*x[5]*x[6] + x[2]*x[3]*x[5] + x[3]*x[4]*x[6]
    G168_RC = [Permutation(6), Permutation(2,4,5), Permutation(2,5,4), Permutation(2,3)(4,5), Permutation(2,4,3), Permutation(2,5,3),
               Permutation(3,4,5), Permutation(2,3,4), Permutation(2,5)(3,4), Permutation(3,5,4), Permutation(2,4)(3,5),
               Permutation(2,3,5), Permutation(3,6)(4,5), Permutation(2,4)(3,6), Permutation(2,5)(3,6), Permutation(1,3,2,6,4,5),
               Permutation(1,3,2,5,6,4), Permutation(1,3,2)(4,6), Permutation(1,3,6,4), Permutation(1,3,6,4,2,5), Permutation(1,3,6,4,5,2),
               Permutation(1,3,5,2,6,4), Permutation(1,3,5)(4,6), Permutation(1,3,5,6,4,2), Permutation(1,3)(2,6,4),
               Permutation(1,3)(2,5)(4,6), Permutation(1,3)(4,5,6), Permutation(1,3,4)(2,6), Permutation(1,3,4,6,2,5),
               Permutation(1,3,4,5,6,2)]
    fix_perm = compute_resolvent_roots(f_roots, F_G168, G168_RC)

    if fix_perm is not None:

        f_roots = permute_list(f_roots, fix_perm)

        F_G21 = lambda x: x[0] * x[1] * x[3] + x[0] * x[1] * x[5] + x[0] * x[2] * x[3] + x[0] * x[2] * x[6] + x[0] * x[4] * x[5] + x[0] * x[4] * x[6] + x[1] * x[2] * x[4] + x[1] * x[2] * x[6] + x[1] * x[3] * x[4] + x[1] * x[5] * x[6] + x[2] * x[3] * x[5] + x[2] * x[4] * x[5] + x[3] * x[4] * x[6] + x[3] * x[5] * x[6]
        G21_RC = [Permutation(6), Permutation(2,6)(4,5), Permutation(1,2)(6,3), Permutation(1,2,3,6)(4,5), Permutation(1,3)(4,5), Permutation(1,3)(2,6), Permutation(1,6,3,2)(4,5), Permutation(1,6)(2,3)]
        fix_perm = compute_resolvent_roots(f_roots, F_G21, G21_RC)
        if fix_perm is not None:

            f_roots = permute_list(f_roots, fix_perm)

            F_G7 = lambda x: x[0]*x[1] + x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6] + x[6]*x[0]
            G7_RC = [Permutation(6), Permutation(1,2,4)(3,6,5), Permutation(1,4,2)(3,5,6)]
            fix_perm = (f_roots, F_G7, G7_RC)
            if fix_perm is not None:
                return 'G7'
            else:
                return 'G21'

        else:
            return 'G168 (PSL_3(2))'

    F_G42 = lambda x: x[0]*x[1]*x[3] + x[0]*x[1]*x[5] + x[0]*x[2]*x[3] + x[0]*x[2]*x[6] + x[0]*x[4]*x[5] + x[0]*x[4]*x[6] + x[1]*x[2]*x[4] + x[1]*x[2]*x[6] + x[1]*x[3]*x[4] + x[1]*x[5]*x[6] + x[2]*x[3]*x[5] + x[2]*x[4]*x[5] + x[3]*x[4]*x[6] + x[3]*x[5]*x[6]
    G42_RC = [Permutation(6), Permutation(2, 6)(4, 5), Permutation(1, 2)(3, 6), Permutation(1, 2, 3, 6)(4, 5), Permutation(6)(1, 3)(4, 5), Permutation(1, 3)(2, 6), Permutation(1, 6, 3, 2)(4, 5),
            Permutation(1, 6)(2, 3), Permutation(6)(2, 4, 5), Permutation(2, 5, 6), Permutation(1, 2, 4, 5)(3, 6), Permutation(1, 2, 5, 3, 6), Permutation(1, 3)(2, 5), Permutation(1, 3)(2, 4, 5, 6),
            Permutation(1, 6, 3, 2, 5), Permutation(1, 6)(2, 4, 5, 3), Permutation(6)(2, 5, 4), Permutation(2, 4, 6), Permutation(1, 2, 5, 4)(3, 6), Permutation(1, 2, 4, 3, 6), Permutation(5)(1, 3)(2, 4),
            Permutation(1, 3)(2, 5, 4, 6), Permutation(1, 6, 3, 2, 4), Permutation(1, 6)(2, 5, 4, 3), Permutation(6)(2, 3)(4, 5), Permutation(2, 3, 6), Permutation(1, 2, 6, 3)(4, 5), Permutation(1, 2, 6),
            Permutation(5)(1, 3, 2), Permutation(1, 3, 6, 2)(4, 5), Permutation(1, 6, 3), Permutation(1, 6)(4, 5), Permutation(6)(2, 4, 3), Permutation(2, 5, 4, 3, 6), Permutation(1, 2, 4, 6, 3),
            Permutation(1, 2, 5, 4, 6), Permutation(1, 3, 2, 5, 4), Permutation(1, 3, 6, 2, 4), Permutation(1, 6, 3)(2, 5, 4), Permutation(1, 6)(2, 4), Permutation(6)(2, 5, 3), Permutation(2, 4, 5, 3, 6),
            Permutation(1, 2, 5, 6, 3), Permutation(1, 2, 4, 5, 6), Permutation(1, 3, 2, 4, 5), Permutation(1, 3, 6, 2, 5), Permutation(1, 6, 3)(2, 4, 5), Permutation(1, 6)(2, 5), Permutation(6)(3, 4, 5),
            Permutation(2, 6)(3, 5), Permutation(1, 2)(3, 4, 5, 6), Permutation(1, 2, 3, 5, 6), Permutation(1, 3, 5), Permutation(1, 3, 4, 5)(2, 6), Permutation(1, 6, 3, 5, 2), Permutation(1, 6)(2, 3, 4, 5),
            Permutation(6)(2, 3, 4), Permutation(2, 3, 5, 4, 6), Permutation(1, 2, 6, 3, 4), Permutation(1, 2, 6)(3, 5, 4), Permutation(1, 3, 5, 4, 2), Permutation(1, 3, 4, 6, 2), Permutation(1, 6, 3, 5, 4),
            Permutation(1, 6)(3, 4), Permutation(6)(2, 5)(3, 4), Permutation(2, 4, 3, 5, 6), Permutation(1, 2, 5)(3, 4, 6), Permutation(1, 2, 4, 6)(3, 5), Permutation(1, 3, 5, 2, 4), Permutation(1, 3, 4)(2, 5, 6),
            Permutation(1, 6, 3, 5)(2, 4), Permutation(1, 6)(2, 5, 3, 4), Permutation(6)(3, 5, 4), Permutation(2, 6)(3, 4), Permutation(1, 2)(3, 5, 4, 6), Permutation(1, 2, 3, 4, 6), Permutation(5)(1, 3, 4),
            Permutation(1, 3, 5, 4)(2, 6), Permutation(1, 6, 3, 4, 2), Permutation(1, 6)(2, 3, 5, 4), Permutation(6)(2, 4)(3, 5), Permutation(2, 5, 3, 4, 6), Permutation(1, 2, 4)(3, 5, 6), Permutation(1, 2, 5, 6)(3, 4),
            Permutation(1, 3, 4, 2, 5), Permutation(1, 3, 5)(2, 4, 6), Permutation(1, 6, 3, 4)(2, 5), Permutation(1, 6)(2, 4, 3, 5), Permutation(6)(2, 3, 5), Permutation(2, 3, 4, 5, 6), Permutation(1, 2, 6, 3, 5),
            Permutation(1, 2, 6)(3, 4, 5), Permutation(1, 3, 4, 5, 2), Permutation(1, 3, 5, 6, 2), Permutation(1, 6, 3, 4, 5), Permutation(1, 6)(3, 5), Permutation(3, 6)(4, 5), Permutation(2, 6, 3),
            Permutation(6)(1, 2)(4, 5), Permutation(6)(1, 2, 3), Permutation(1, 3, 6), Permutation(1, 3, 2, 6)(4, 5), Permutation(1, 6, 2), Permutation(1, 6, 2, 3)(4, 5), Permutation(2, 4)(3, 6),
            Permutation(2, 5, 4, 6, 3), Permutation(6)(1, 2, 4), Permutation(6)(1, 2, 5, 4, 3), Permutation(1, 3, 6)(2, 5, 4), Permutation(1, 3, 2, 4, 6), Permutation(1, 6, 2, 5, 4), Permutation(1, 6, 2, 4, 3),
            Permutation(2, 5)(3, 6), Permutation(2, 4, 5, 6, 3), Permutation(6)(1, 2, 5), Permutation(6)(1, 2, 4, 5, 3), Permutation(1, 3, 6)(2, 4, 5), Permutation(1, 3, 2, 5, 6), Permutation(1, 6, 2, 4, 5),
            Permutation(1, 6, 2, 5, 3)]
    fix_perm = compute_resolvent_roots(f_roots, F_G42, G42_RC)
    if fix_perm is not None:

        f_roots = permute_list(f_roots, fix_perm)

        F_G14 = lambda x: x[0]*x[1] + x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[5] + x[5]*x[6] + x[6]*x[0]
        G14_RC = [Permutation(7), Permutation(1,2,4)(3,6,5), Permutation(1,4,2)(3,5,6)]
        fix_perm = compute_resolvent_roots(f_roots, F_G14, G14_RC)
        if fix_perm is not None:
            return 'G14'
        else:
            return 'G42'

    if is_square(discriminant(f)):
        return 'S7'
    else:
        return 'A7'

