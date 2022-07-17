# 6_3_PROB
import scripts.Rosalind as rs
from math import log10

### read data ###
with open("../../inputs/stronghold/rosalind_prob.txt", "r") as f:
    seq = f.readline().rstrip()
    arr_A = [float(x) for x in f.readline().rstrip().split()]

seq = rs.nucleotides(seq)
counts = seq.count()

### calculate array B ###
at = counts["A"] + counts["T"]
gc = counts["G"] + counts["C"]

arr_B = []
for gc_p in arr_A:
    at_p = 1 - gc_p
    b_prob = log10((at_p / 2) ** at)
    b_prob += log10((gc_p / 2) ** gc)
    arr_B.append(b_prob)

### print ###
print(" ".join(map(str, arr_B)))
