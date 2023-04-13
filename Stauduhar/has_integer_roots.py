from sympy import symbols, roots

def has_integer_roots(f):
    x = symbols('x')
    f_roots = roots(f, x)

    for root, multiplicity in f_roots.items():
        if root.is_integer:
            return True
    return False
