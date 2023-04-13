from sympy import symbols, nroots, Poly
from permute_list import permute_list

def compute_resolvent_roots(f, F, perm_lst, tol=10e-10):

    x = symbols('x')

    f_roots = nroots(f, n=15)
    R_int_roots = set()
    integer_root, integer_perm = None, None

    for perm in perm_lst:

        permuted_roots = permute_list(nroots, perm)
        root_of_F = F(permuted_roots)

        if abs(root_of_F.round() - root_of_F) < tol:
            if root_of_F.round() in R_int_roots:
                raise ValueError('The resolvent has a repeated integer root')
            else:
                R_int_roots.add(root_of_F.round())
                integer_root = root_of_F.round()
                integer_perm = perm

    return integer_perm

