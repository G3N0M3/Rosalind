# 5_2_FIBD
"""
See note
"""
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_fibd.txt", "r") as f:
    n, k = map(int, f.readline().rstrip().split())

res = rs.fib_death(gen=n, lifespan=k, litter=1, init=1)
print(sum(res))
