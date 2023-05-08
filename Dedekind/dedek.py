from sympy import prime, Poly, symbols, factor_list, degree, discriminant
from collections import Counter

class Partition:
    def __init__(self, decomp, nprimes):
        self.decomp = decomp
        self.nprimes = nprimes
    def addprime(self, n=1):
        self.nprimes += n
    def equal(self, deco):
        return self.decomp == deco

def induced_partition(f, p):

    return Counter([degree(factor[0]) for factor in factor_list(f, modulus=p)[1]])

def add_partition_to_list(lst, partition):

    for part in lst:
        if part.equal(partition):
            part.addprime()
            return

    lst.append(Partition(partition, 1))
    return
def discover_gal(f, max_n_p):

    n_prime = 1
    total_primes = 0
    disc = discriminant(f)
    partition_list = []

    while(n_prime <= max_n_p):

        p = prime(n_prime)

        if(disc % p == 0):
            n_prime += 1
            continue
        else:
            total_primes += 1

        partition = induced_partition(f, p)
        add_partition_to_list(partition_list, partition)

        n_prime += 1

    return partition_list, total_primes

def counter_to_list(count):
    count = dict(count)
    lst = []
    for key, val in count.items():
        lst += [key]*val
    return lst




