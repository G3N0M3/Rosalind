# 7_12_REAR
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_rear.txt", "r") as f:
    permutation_set, _append = [], []
    for _line in f:
        line = _line.rstrip()
        if line == "":
            permutation_set.append(_append)
            _append = []
        else:
            _append.append(line.split())
    permutation_set.append(_append)

distances = []
for item in permutation_set:
    print(f"Comparing {item[0]} and {item[1]}")
    dist = rs.reversal_distance(item[0], item[1])
    distances.append(dist)

print(" ".join(str(x) for x in distances))
