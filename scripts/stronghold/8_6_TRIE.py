# 8_6_TRIE
import scripts.Rosalind as rs
import pandas as pd

### read data ###
with open("../../inputs/stronghold/rosalind_trie.txt", "r") as f:
    seqs = []
    for line in f:
        seq = line.rstrip()
        seqs.append(seq)

print(seqs)
division = rs.divide_by_first(seqs)
print(division)
