# 9_4_DBRU

### read data ###
with open("../../inputs/stronghold/rosalind_dbru.txt", "r") as f:
    k_plus_mers = {seq.rstrip() for seq in f.readlines()}

print(k_plus_mers)
