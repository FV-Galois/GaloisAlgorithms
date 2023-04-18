from degs import *
from sympy import symbols, Poly

x = symbols('x')
f = Poly(x**4 + 8*x + 12, x)

print(deg4(f))