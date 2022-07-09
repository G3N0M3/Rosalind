# 5_4_LIA
"""
See note
"""
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_lia.txt", "r") as f:
    k, N = map(int, f.readline().rstrip().split())

res = rs.binomial_sum_limit(p=.25, n=2**k, limit=N-1, coef=False)
print(1 - res)
