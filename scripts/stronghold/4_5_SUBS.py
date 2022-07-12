# 4_5_SUBS
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_subs.txt", "r") as f:
    sup, sub = map(lambda x: x.rstrip(), list(f))

### substring search ###
sub_idxs = rs.substring_idx(sup, sub)
for idx in sub_idxs:
    print(idx, end=" ")
