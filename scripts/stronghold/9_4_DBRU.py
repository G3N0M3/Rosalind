# 9_4_DBRU
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_dbru.txt", "r") as f:
    seq = f.readline().rstrip()
    k_plus_mers = {seq.rstrip() for seq in f.readlines()}
    k_plus_mers = k_plus_mers | {seq}
    k = len(seq) - 1

### get reverse complement and merge ###
k_plus_mers_reverse = {rs.reverse_complement(seq) for seq in k_plus_mers}
k_plus_mers_total = k_plus_mers | k_plus_mers_reverse

### get output list ###
res = [(r[:k], r[1:k + 1]) for r in k_plus_mers_total]

### print ###
for item in res:
    print(f"({item[0]}, {item[1]})")
