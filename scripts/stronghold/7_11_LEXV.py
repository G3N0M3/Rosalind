# 7_11_LEXV
import scripts.Rosalind as rs
from itertools import product

### read data ###
with open("../../inputs/stronghold/rosalind_lexv.txt", "r") as f:
    alpha = list(f.readline().rstrip().split())
    n = int(f.readline().rstrip())

### building criteria dictionary ###
lex_crit = {}
for idx in range(len(alpha)):
    lex_crit[alpha[idx]] = idx

### Every permutation of length 1 to n ###
prods = []
for _len in range(1, n + 1):
    words = product(alpha, repeat=_len)
    for word in words:
        _word = "".join(word)
        prods.append(_word)

### Bubble sort using defined order method ###
prod_idx = len(prods)
for idx1 in range(prod_idx - 1):
    print(f"Working for index ")
    for idx2 in range(idx1 + 1, prod_idx):
        prods[idx1], prods[idx2] = rs.order_lex(prods[idx1], prods[idx2], criteria=lex_crit)

### print ###
for item in prods:
    print(item)
