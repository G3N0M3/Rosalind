# 5_9_ORF
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_orf.txt", "r") as f:
    name = f.readline().rstrip()
    seq = "".join(x.rstrip() for x in f.readlines())
    seq = rs.nucleotides(seq)

orf_li = seq.orf(rev_comp=True)

for orf in orf_li:
    print(orf)
