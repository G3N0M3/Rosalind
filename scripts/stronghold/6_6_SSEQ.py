# 6_6_SSEQ
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_sseq.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))

seqs = list(fastd.values())
sup = seqs[0]
sub = seqs[1]

res = rs.spliced_motif(sup, sub, case=1, base=1)

print(" ".join(map(str, res[0])))
