# 8_9_EDIT
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_edit.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seqs = list(fastd.values())

print(seqs)
