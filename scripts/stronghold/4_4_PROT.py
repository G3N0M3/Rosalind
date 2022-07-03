# 4_4_PROT
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_prot.txt", "r") as f:
    seq = rs.read_seq(f)

    aa_seq = rs.nucleotides(seq).translate()

print(aa_seq)
