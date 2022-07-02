# 3_1_REVC
import scripts.Rosalind as rs

with open("../../inputs/stronghold/rosalind_revc.txt", "r") as f:
    seq = rs.read_seq(f)
    seq = rs.nucleotides(seq)
    seq_comp = seq.complement(reverse=True)
    print(seq_comp)
