from degs import *
from sympy import symbols, Poly

x = symbols('x')
f = Poly(x**4 + x**3 + x**2 + x + 1, x)

print(deg_4(f))