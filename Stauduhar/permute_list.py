from sympy.combinatorics import Permutation

def permute_list(lst, perm):

    permuted_lst = [None] * len(lst)
    for i in range(len(lst)):
        try:
            permuted_lst[perm(i)] = lst[i]
        except:
            permuted_lst[i] = lst[i]
    return permuted_lst

