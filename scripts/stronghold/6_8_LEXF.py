# 6_8_LEXF
from itertools import product

### read data ###
with open("../../inputs/stronghold/rosalind_lexf.txt", "r") as f:
    symbols = f.readline().split()
    n = int(f.readline().rstrip())

res = []
for item in product(symbols, repeat=n):
    res.append("".join(item))

for item in sorted(res):
    print(item)
