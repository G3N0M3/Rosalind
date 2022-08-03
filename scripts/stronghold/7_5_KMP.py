# 7_5_KMP
import scripts.Rosalind as rs


### read data ###
with open("../../inputs/stronghold/rosalind_kmp.txt", "r") as f:
    seq_id = f.readline().rstrip()
    seq = rs.read_seq(f.readlines())

failure_array = [-1, 0]
"""
array[0] initial process for 1-base numbering
array[1] = 0
"""
### complete failure array ###
"""
Theoretically works properly, but takes too much time
"""
for idx in range(1, len(seq)):
    print(f"Checking index {idx + 1} / {len(seq)}")
    check_seq = seq[:idx + 1]
    match_length = rs.longest_prefix_suffix_match(check_seq)
    failure_array.append(match_length)

### print ###
print(" ".join([str(x) for x in failure_array[1:]]))
