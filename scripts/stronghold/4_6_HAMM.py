# 4_6_HAMM
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_hamm.txt", "r") as f:
    s, t = map(lambda x: x.rstrip(), list(f))

### calculate Hamming distance
hamm = rs.hamm_distance(s, t)

print(hamm)
