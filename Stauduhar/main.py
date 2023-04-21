from degs import *
from sympy import symbols, Poly

x = symbols('x')
f = Poly(x**7 - 7*x**3 + 14*x**2 - 7*x + 1, x)

print(deg_7(f))