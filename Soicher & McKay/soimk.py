from sympy import *
from random import randint
from itertools import combinations
from random import randint

def round(x) -> int:
	'''
	sign = int(x/abs(x))
	x = abs(x)
	floor = int(x)
	return sign * floor if x < floor + 0.5 else sign * (floor + 1)
	'''
	return re(x).round()

def compute_resolvent(F, f, r):
	x = Symbol('x')
	f_roots = nroots(f, n=15)

	perms = combinations(f_roots, r)	#Computes the action of S_n on the set of roots of f

	R_roots = []
	for p in perms:	#Computes the F_i, ie, the roots of R(F,f)
		R_roots.append(F(p))

	R = Poly(1, x)
	for root in R_roots:	#Computes R(F,f)
		R *= Poly(x - root, x, domain='CC')

	coeffs = R.coeffs()
	R = Poly(0, x)	#Rounds the coefficients of R(F,f)
	for i, coef in enumerate(reversed(coeffs)):
		R += Poly(round(coef)*x**i, x)

	return R

def tsc_trans(f):

	def tsc_trans_aux(f): #Computes a possible Tschirnhaus transformation

		x = Symbol('x')
		y = Symbol('y')

		f = f.subs(x, y)

		n = degree(f)
		deg = randint(1, n - 1)

		g = Poly(1, y)

		for i in range(deg): #generates a random polynomial
			g *= Poly(y - randint(-10, 10))

		u = resultant(f, Poly(x, x) - g, gens=y)

		u = u.subs(y, x)
		udif = diff(u, x)

		d = gcdex(u, udif, gens=x)[2] #computes gcd(u, u')

		return u, degree(d)

	u, dd = tsc_trans_aux(f)

	while dd != 0: #Generate a new transformation until gcd(u, u') is constant
		u, d = tsc_trans_aux(f)

	return u

def has_repeated_roots(f):
	roots = nroots(f, n=15, maxsteps=10000)

	for i in range(len(roots)):
		for j in range(i+1, len(roots)):
			if abs(roots[i] - roots[j]) < 10e-8:
				return True

	return False

def is_irred(f):
	return f.factor_list()[0] == 1

def generate_a():
	a = []
	while(len(a) != 2):
		b = randint(-100, 100)
		if b != 0 and b not in a:
			a.append(b)
	return a

def orbit_length_partition_set(F, f, r):

	R = compute_resolvent(F, f, r)

	while(has_repeated_roots(R)):
		f = tsc_trans(f)
		R = compute_resolvent(F, f, r)

	facts = [fact[0] for fact in factor_list(R)[1]]
	OLP = [degree(fac) for fac in facts]

	return OLP

def orbit_length_partition_seq(f, r):

	a = generate_a()
	F = seq_funct(a)
	R = compute_resolvent(F, f, r)

	while(has_repeated_roots(R)):
		a = generate_a()
		F = seq_funct(a)
		R = compute_resolvent(F, f, r)

	facts = [fact[0] for fact in factor_list(R)[1]]
	OLP = [degree(fac) for fac in facts]

	return OLP

def L2(x):
	return x[0] + x[1]

def L3(x):
	return x[0] + x[1] + x[2]

def seq_funct(a):
	def L2(x):
		return a[0]*x[0] + a[1]*x[1]
	return L2

x = Symbol('x')

#f = Poly(x**7 - 14*x**5 + 56*x**3 - 56*x + 22, x)
f = Poly(x**6 + 108, x)

op2s = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)
op3s = orbit_length_partition_set(L3, f, 3)
op2sq = orbit_length_partition_seq(f, 2)

print(op2s)
print(op3s)
print(op2sq)