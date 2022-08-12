# 8_7_SCSP
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_scsp.txt", "r") as f:
    seq1 = f.readline().rstrip()
    seq2 = f.readline().rstrip()

scs = rs.shortest_common_supersequence(seq1, seq2)

print(scs)
