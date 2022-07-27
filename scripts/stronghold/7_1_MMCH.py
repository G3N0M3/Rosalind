# 7_1_SIGN
import scripts.Rosalind as rs
from math import factorial

### read data ###
with open("../../inputs/stronghold/rosalind_mmch.txt", "r") as f:
    seq_id = f.readline().rstrip()
    seq = rs.read_seq(f.readlines())
    seq = rs.nucleotides(seq)

### calculate ###
counts = seq.count()
t, T = sorted([counts["A"], counts["U"]])
r, R = sorted([counts["G"], counts["C"]])
res = (factorial(T) // factorial(T - t)) * (factorial(R) // factorial(R - r))
_res = (factorial(T) / factorial(T - t)) * (factorial(R) / factorial(R - r))
"""
See the url below for the reason for using // instead of /
https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
When T, t, R, r = 21, 21, 27, 26
// -> 556322599406617581545815578704263249920000000000
 / -> 556322599406617603788842354602141456698271858688

Can also use math.perm instead of factorial calculation
"""

print(res)
