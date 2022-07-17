# 6_2_PPER
from math import factorial

### read data ###
with open("../../inputs/stronghold/rosalind_pper.txt", "r") as f:
    n, k = [int(x) for x in f.readline().split()]

### calculate partial permutation ###
"""
nPk = n! / (n-k)!
However, overflow while the calculation of factorial might give an overflow problem
Thus, will not use the factorial function from the math library
"""
res = 1
for num in range((n - k + 1), n + 1):
    res *= num
    res %= 1000000

### print
print(res)
