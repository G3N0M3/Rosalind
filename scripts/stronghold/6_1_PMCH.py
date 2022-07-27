# 6_1_PMCH
import scripts.Rosalind as rs
from math import factorial

### read data ###
with open("../../inputs/stronghold/rosalind_pmch.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seq = list(fastd.values())[0]
    seq = rs.nucleotides(seq)

"""
See note
"""
counts = seq.count()
n_au = counts["A"] + counts["U"]  # number of A and U
n_cg = counts["C"] + counts["G"]  # number of C and G
# since the number of A and U, or C and G is equal,
# dividing by two is not necessary, but fully written for understandings
res = factorial(int(n_au / 2)) * factorial(int(n_cg / 2))

print(res)
