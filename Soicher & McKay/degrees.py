from orbit_length_functions import orbit_length_partition_set, orbit_length_partition_seq
from sympy import *

def deg2(f):
	return 'S2'

def deg3(f):
	x = Symbol('x')
	disc = discriminant(f, gens=x)
	return '+A3' if sqrt(disc).is_rational else 'S3'

def deg4(f):
	x = Symbol('x')

	olp2 = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)

	if olp2 == (2, 2, 2): return '+V4'
	if olp2 == (6,):
		disc = discriminant(f, gens=x)
		return '+A4' if sqrt(disc).is_rational else 'S4'

	olp2seq = () #TODO

	if olp2seq == (4, 4, 4): return 'Z4'
	return 'D4'

def deg5(f):
	x = Symbol('x')

	olp2 = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)

	if olp2 == (10,):
		disc = discriminant(f, gens=x)
		if sqrt(disc).is_rational: return '+A5'

		#TODO falta diferenciar entre F20 y S5

	olp2seq = () #TODO

	if olp2seq == (5,5,5,5): return '+Z5'
	return '+D5'

def deg6(f):
	#TODO
	return 0

def deg7(f):
	x = Symbol('x')

	olp2 = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)
	olp3 = orbit_length_partition_set(lambda x: x[0] + x[1] + x[2], f, 3)

	if olp2 == (7,7,7):
		if olp3 == (7,7,7,7,7): return '+Z7'
		else: return 'D7'

	if olp3 == (7,7,21): return '+F21'
	if olp3 == (14, 21): return 'F42'
	if olp3 == (7,28): return '+PSL3(2)'

	disc = discriminant(f, gens=x)
	print(olp2 == (21,))
	return '+A7' if sqrt(disc).is_rational else 'S7'