from sympy import *
from itertools import combinations
from random import randint
from collections import Counter

def orbit_length_partition_set(F, f, r):

	R = compute_resolvent(F, f, r, False)

	while(has_repeated_roots(R)):
		ft = tsc_trans(f)
		R = compute_resolvent(F, ft, r, False)

	facts = [fact[0] for fact in factor_list(R)[1]]
	OLP = Counter([degree(fac) for fac in facts])

	return OLP

def orbit_length_partition_seq(f, r):

	a = generate_a(r)
	F = seq_funct(a, r)
	R = compute_resolvent(F, f, r, True)

	while(has_repeated_roots(R)):
		a = generate_a(r)
		F = seq_funct(a, r)
		R = compute_resolvent(F, f, r, True)

	facts = [fact[0] for fact in factor_list(R)[1]]
	OLP = Counter([degree(fac) for fac in facts])

	return OLP

def compute_resolvent(F, f, r, ordered):
	x = Symbol('x')
	f_roots = nroots(f, n=15)

	perms = combinations(f_roots, r)	#Computes the action of S_n on the set of roots of f

	R_roots = []
	for p in perms:	#Computes the F_i, ie, the roots of R(F,f)
		R_roots.append(F(p))
		if ordered: R_roots.append(F(p[::-1]))

	R = Poly(1, x)
	for root in R_roots:	#Computes R(F,f)
		R *= Poly(x - root, x, domain='CC')

	coeffs = R.all_coeffs()
	R = Poly(0, x)	#Rounds the coefficients of R(F,f)
	for i, coef in enumerate(reversed(coeffs)):
		R += Poly(round(coef)*x**i, x)

	return R

def tsc_trans(f):
	raise NotImplementedError("The tsc trans does not work")

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
			if abs(roots[i] - roots[j]) < 10e-10:
				return True

	return False

def is_irred(f):
	return f.factor_list()[0] == 1

def generate_a(n):

	a = set()

	while(len(a) < n):
		a.add(randint(1,5))

	return list(a)
def round(x) -> int:
	return re(x).round()
def seq_funct(a, r):
	return lambda x: a[0]*x[0] + a[1]*x[1] if r == 2 else lambda x: a[0]*x[0] + a[1]*x[1] + a[2]*x[2]


