# 7_3_SSET

### read data ###
with open("../../inputs/stronghold/rosalind_sset.txt", "r") as f:
    n = int(f.readline().rstrip())

### calculate ###
"""
Since there is two choices for each element (included, excluded)
to create a subset, the number of total subsets for a set is
2^n
"""
# modulo was applied here to avoid calculating large numbers
res = 1
for i in range(n):
    res *= 2
    res %= 1000000

print(res)
