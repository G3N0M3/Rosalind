# 7_6_PDST
import scripts.Rosalind as rs
import pandas as pd


### read data ###
with open("../../inputs/stronghold/rosalind_pdst.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seqs = list(fastd.values())

distance_matrix = rs.p_distance_matrix(seqs)

### print ###
for idx in range(1, len(seqs) + 1):
    row = distance_matrix.loc[:, idx]
    print(" ".join([str(x) for x in row]))
