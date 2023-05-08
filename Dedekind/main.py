from sympy import symbols, Poly
from dedek import discover_gal, Partition, counter_to_list

x = symbols('x')
f = Poly(x**7 - 7*x**3 + 14*x**2 - 7*x + 1, x)

cycle_types, n = discover_gal(f, 200)
for cycle_type in cycle_types:
    print(f"Cycle Type: {counter_to_list(cycle_type.decomp)} Estimated density: {cycle_type.nprimes / n}")