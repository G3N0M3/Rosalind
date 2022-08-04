# 7_8_kmer
import scripts.Rosalind as rs
from itertools import product

### read data ###
with open("../../inputs/stronghold/rosalind_kmer.txt", "r") as f:
    seq_id = f.readline().rstrip()
    seq = rs.read_seq(f.readlines())

# configuring "k"
k = 4

# generating kmer dictionary
kmers = sorted(["".join(x) for x in product("ATGC", repeat=k)])
kmer_d = {x: 0 for x in kmers}

# counting kmer appearance in sequence
for start_idx in range(len(seq) - k + 1):
    kmer = seq[start_idx:start_idx + 4]
    kmer_d[kmer] += 1

### print ###
print(" ".join([str(x) for x in kmer_d.values()]))
