# 8_3_SETO
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_scsp.txt", "r") as f:
    seq1 = f.readline().rstrip()
    seq2 = f.readline().rstrip()

scs = rs.shortest_common_supersequence(seq1, seq2)
"""
Function currently generates wrong output
See reference in function docstring
"""

print(scs)
