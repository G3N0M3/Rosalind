# 1_1_DNA
import scripts.Rosalind as rs

with open("../../inputs/stronghold/rosalind_dna.txt", "r") as f:
    seq = rs.read_seq(f)
    seq = rs.nucleotides(seq)
    counts = seq.count()

    for base in "ACGT":
        print(counts[base], end=" ")
