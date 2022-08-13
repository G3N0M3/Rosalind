# 9_3_INDC
from math import comb
from math import log10

### read data ###
with open("../../inputs/stronghold/rosalind_indc.txt", "r") as f:
    n = int(f.readline().rstrip())

res = []
_sum = 0
for k in range(2 * n, 0, -1):
    _sum += comb(2 * n, k) * (.5 ** k) * (.5 ** (2 * n - k))
    res.append(log10(_sum))

for cal in res[::-1]:
    print(f"{cal:.4f}", end=" ")
