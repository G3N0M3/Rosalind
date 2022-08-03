# 7_7_rstr
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_rstr.txt", "r") as f:
    given = f.readline().rstrip().split()
    n, gc = int(given[0]), float(given[1])
    seq = rs.read_seq(f.readlines())

### calculate probability ###
# create probability dictionary for each nucleotide
_g, _c = gc / 2, gc / 2
_a, _t = (1 - gc) / 2, (1 - gc) / 2
nt_prob = {"A": _a, "T": _t, "G": _g, "C": _c}
# calculate motif probability
motif_prob = 1
for nt in seq:
    motif_prob *= nt_prob[nt]
# calculate answer
ans = 1 - (1 - motif_prob) ** n

print(ans)
