from sympy import *
from permute_list import permute_list
from random import randint

def compute_resolvent_roots(f_roots, F, perm_lst, tol=10e-10): #Non exact approach

    x = symbols('x')

    R_int_roots = set()
    integer_perm = None

    for perm in perm_lst:

        permuted_roots = permute_list(f_roots, perm)
        root_of_F = N(F(permuted_roots),20)

        if abs(im(root_of_F)) < tol and abs(re(root_of_F).round() - re(root_of_F)) < tol:
            int_root = re(root_of_F).round()
            if int_root in R_int_roots:
                raise ValueError('The resolvent has a repeated integer root')
            else:
                R_int_roots.add(int_root)
                integer_perm = perm

    return integer_perm

def return_integer_permutation(f_roots, F, perm_lst): #Exact approach

    roots = [F(permute_list(f_roots, perm)) for perm in perm_lst]
    res_poly, a_0 = generate_resolvent(roots)

    if(has_repeated_roots(res_poly)):
        raise ValueError("Has repeated roots")

    int_root = has_integer_root(res_poly, a_0)
    print(res_poly)

    if int_root is None:
        return None

    for i in range(len(perm_lst)):

        perm_lst_partial = perm_lst[:i] + perm_lst[i+1:]
        roots_partial = [F(permute_list(f_roots, perm)) for perm in perm_lst_partial]
        res_poly_partial, _ = generate_resolvent(roots_partial)

        if res_poly_partial.eval(int_root) != 0:
            return perm_lst[i]

    raise ArithmeticError("Something Failed")

def has_repeated_roots(f):
    roots = nroots(f, n=15, maxsteps=10000)

    for i in range(len(roots)):
        for j in range(i+1, len(roots)):
            if abs(roots[i] - roots[j]) < 10e-10:
                return True

    return False

def has_integer_root(res_poly, coef): #coef is the coeficient a_0

    candidates = divisors(coef) if coef != 0 else [0]
    roots = set()

    for candidate in candidates:
        if res_poly.eval(candidate) == 0:
            if candidate in roots:
                raise ValueError('The resolvent has a repeated integer root')
            else:
                roots.add(candidate)

    return roots.pop() if roots else None

def generate_resolvent(roots):

    x = symbols('x')
    res_poly = Poly(x**len(roots), x)

    coefs = coefs_from_roots(roots)
    for i, coef in enumerate(reversed(coefs[1:])):
        res_poly += Poly(coef * x ** i, x)

    return res_poly, coefs[-1]

def coefs_from_roots(roots):

    return [1] + [re(N(coef, 10)).round() for coef in rec_coefs_from_roots(roots, [1])[1:]]
def rec_coefs_from_roots(roots, coefs): #coefs are in descending order a_n,a_n-1, ..., a_0

    if len(roots) == 0:
        return coefs

    root = roots[0]
    coefs = [coefs[0]] + coefs

    for i in range(1, len(coefs)-1):
        coefs[i] = coefs[i+1] - root*coefs[i]
    coefs[-1] = -root*coefs[-1]

    return rec_coefs_from_roots(roots[1:], coefs)
