from sympy import symbols, nroots, Poly, re, im
from permute_list import permute_list

def compute_resolvent_roots(f_roots, F, perm_lst, tol=10e-10):

    x = symbols('x')

    R_int_roots = set()
    integer_root, integer_perm = None, None

    for perm in perm_lst:

        permuted_roots = permute_list(f_roots, perm)
        root_of_F = F(permuted_roots)

        if abs(im(root_of_F)) < tol and abs(re(root_of_F).round() - root_of_F) < tol:
            int_root = re(root_of_F).round()
            if int_root in R_int_roots:
                raise ValueError('The resolvent has a repeated integer root')
            else:
                R_int_roots.add(int_root)
                integer_root = int_root
                integer_perm = perm

    return integer_perm

