# 5_12_PERM
from itertools import permutations as perm

### read data ###
with open("../../inputs/stronghold/rosalind_perm.txt", "r") as f:
    n = int(f.readline())

print(n)
for item in perm(range(1, n+1)):
    res = " ".join([str(x) for x in item])
    print(res)
