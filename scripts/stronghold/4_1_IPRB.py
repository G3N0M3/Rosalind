# 4_1_IPRB
from math import comb

### read data ###
with open("../../inputs/stronghold/rosalind_iprb.txt", "r") as f:
    k, m, n = map(int, f.readline().rstrip().split())
    total = k + m + n

    ## calculation
    res = 0
    # AA x AA
    res += comb(k, 2) / comb(total, 2) * 1
    # AA x Aa
    res += k * m / comb(total, 2) * 1
    # AA x aa
    res += k * n / comb(total, 2) * 1
    # Aa x Aa
    res += comb(m, 2) / comb(total, 2) * .75
    # Aa x aa
    res += m * n / comb(total, 2) * .5

    print(res)
