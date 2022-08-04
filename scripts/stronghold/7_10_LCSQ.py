# 7_10_LCSQ
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_lcsq.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seqs = list(fastd.values())

lcs = rs.longest_common_subsequence(seqs[0], seqs[1])
"""
Function took about 2 minutes to calculate lcs
See reference in function docstring for improvement on calculation
"""

print(lcs)
