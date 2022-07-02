# 4_2_FIB
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_fib.txt", "r") as f:
    n, k = map(int, f.readline().rstrip().split())

    res = rs.fib_rabbits(gen=n, litter=k, init=1)
    print(res)
