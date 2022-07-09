# 5_3_MRNA
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_mrna.txt", "r") as f:
    seq = rs.read_seq(f)
    seq = rs.proteins(seq)
    var = seq.reverse_translate_var()

print(var % 1000000)
