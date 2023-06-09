from degs import *
from sympy import symbols, Poly

x = symbols('x')
f = Poly(x**5 + 2, x)

print(deg_5(f))