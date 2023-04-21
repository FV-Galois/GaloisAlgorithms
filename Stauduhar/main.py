from degs import *
from sympy import symbols, Poly

x = symbols('x')
f = Poly(x**3 + 2, x)

print(deg_3(f))