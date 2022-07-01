# 2_1_RNA
import scripts.Rosalind as rs

with open("../../inputs/stronghold/rosalind_rna.txt", "r") as f:
    seq = rs.read_seq(f)
    seq = rs.nucleotides(seq)
    seq_ts = seq.transcribe()
    print(seq_ts)