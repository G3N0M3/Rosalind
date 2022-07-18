# 6_5_LONG
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_long.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seqs = list(fastd.values())

res = rs.fragment_assembly(seqs)
print(res)
