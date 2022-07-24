# 6_9_LGIS
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_lgis.txt", "r") as f:
    n = int(f.readline().rstrip())
    _pi = [int(x) for x in f.readline().rstrip().split()]

res_asc = rs.longest_pattern_subsequence(_pi, ascending=True)
res_desc = rs.longest_pattern_subsequence(_pi, ascending=False)

print(res_asc)
print(res_desc)
