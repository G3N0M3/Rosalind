# 5_6_GRPH
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_grph.txt", "r") as f:
    fastd = {}
    for name, seq in rs.parse_fasta(f):
        fastd[name[1:]] = seq

res = []
names = list(fastd.keys())
for name1 in names:
    seq1 = fastd[name1]
    for name2 in names:
        if name1 != name2:
            seq2 = fastd[name2]
            if rs.overlap(seq1, seq2, n=3):
                res.append((name1, name2))

for item in res:
    print(item[0], item[1])
