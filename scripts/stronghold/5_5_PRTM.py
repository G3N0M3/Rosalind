# 5_5_PRTM
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_prtm.txt", "r") as f:
    seq = rs.read_seq(f)
    seq = rs.proteins(seq)

print(seq.weight())
