# 5_11_LCSM
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_lcsm.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seqs = list(fastd.values())

common = rs.shared_motif(seqs)
print(common)
