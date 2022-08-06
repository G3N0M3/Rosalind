# 8_2_ASPC
from math import comb

### read data ###
with open("../../inputs/stronghold/rosalind_aspc.txt", "r") as f:
    n, m = map(int, f.readline().rstrip().split())

### calculation ###
"""
Sum(k from m to n) nCk
"""
res = 0
for k in range(m, n + 1):
    res += comb(n, k)
    res %= 1000000

print(res)
