from orbit_length_functions import orbit_length_partition_set, orbit_length_partition_seq
from collections import Counter
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

	if olp2 == Counter((2, 2, 2)): return '+V4'
	if olp2 == Counter((6,)):
		disc = discriminant(f, gens=x)
		return '+A4' if sqrt(disc).is_rational else 'S4'

	olp2seq = orbit_length_partition_seq(f, 2)

	if olp2seq == Counter((4, 4, 4)): return 'Z4'
	return 'D4'

def deg5(f):
	x = Symbol('x')

	olp2 = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)

	if olp2 == Counter((10,)):
		disc = discriminant(f, gens=x)
		if sqrt(disc).is_rational: return '+A5'
		else: return 'F20 or S5'


	olp2seq = orbit_length_partition_seq(f, 2)

	if olp2seq == Counter((5,5,5,5)): return '+Z5'
	return '+D5'

def deg6(f):

	olp2 = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)
	olp3 = orbit_length_partition_set(lambda x: x[0] + x[1] + x[2], f, 3)
	olp2seq = orbit_length_partition_seq(f, 2)

	if olp2 == Counter((3,6,6)):
		if olp3 == Counter((2,6,6,6)): return 'Z6'
		else: return 'D6'

	if olp2 == Counter((3,3,3,6)): return 'S3'

	if olp2 == Counter((3,12)):
		if olp3 == Counter((4,4,6,6)): return 'A4'
		if olp3 == Counter((6,6,8)): return 'G24'
		if olp3 == Counter((4,4,12)): return '+S4/V4'
		return 'G48 or S4/Z4'

	if olp2 == Counter((6,9)):
		if olp2seq == Counter((6,6,18)): return 'G18'
		return 'G^1_36 or G^2_36 or G72'

	if olp3 == Counter((10,10)): return 'PSL2(5)'
	return 'S6 or A6 or PGL2(5)'

def deg7(f):
	x = Symbol('x')

	olp2 = orbit_length_partition_set(lambda x: x[0] + x[1], f, 2)
	olp3 = orbit_length_partition_set(lambda x: x[0] + x[1] + x[2], f, 3)

	if olp2 == Counter((7,7,7)):
		if olp3 == Counter((7,7,7,7,7)): return '+Z7'
		else: return 'D7'

	if olp3 == Counter((7,7,21)): return '+F21'
	if olp3 == Counter((14, 21)): return 'F42'
	if olp3 == Counter((7,28)): return '+PSL3(2)'

	disc = discriminant(f, gens=x)
	return '+A7' if sqrt(disc).is_rational else 'S7'