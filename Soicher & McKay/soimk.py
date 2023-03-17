from sympy import *
from degrees import *

x = Symbol('x')

#f = Poly(x**7 - 14*x**5 + 56*x**3 - 56*x + 22, x)
f = Poly(x**7 - 7*x**3 + 14*x**2 - 7*x + 1, x)

print(deg7(f))
'''
op2s = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)
op3s = orbit_length_partition_set(L3, f, 3)
op2sq = orbit_length_partition_seq(f, 2)

print(op2s)
print(op3s)
print(op2sq)
'''