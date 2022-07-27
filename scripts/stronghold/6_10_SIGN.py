# 6_10_SIGN
from itertools import permutations as perm

### read data ###
with open("../../inputs/stronghold/rosalind_sign.txt", "r") as f:
    n = int(f.readline().rstrip())
num_list = range(1, n + 1)
perm_list = [list(x) for x in perm(num_list, n)]

### make permutation list ###
res = []
for prm in perm_list:
    init_num = prm.pop(0)
    _res = [[init_num], [-init_num]]
    while len(prm) != 0:
        num = prm.pop(0)
        temp = []
        for item in _res:
            temp.append(item + [num])
            temp.append(item + [-num])
            _res = temp
    res += _res

### print ###
print(len(res))
for item in res:
    print(" ".join([str(x) for x in item]))
