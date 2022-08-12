# 8_9_EDIT
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_edit.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    strs = list(fastd.values())

dist = rs.edit_distance(strs[0], strs[1])
"""
The function took about 2 minutes to calculate
input data given at 22.08.12
length of str1 -> 880
length of str2 -> 858
"""

print(dist)
