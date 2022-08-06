# 8_3_SETO

### read data ###
with open("../../inputs/stronghold/rosalind_seto.txt", "r") as f:
    n = int(f.readline().rstrip())
    set_U = set(range(1, n + 1))
    set_A = set(map(int, f.readline().strip("{}\n").split(", ")))
    set_B = set(map(int, f.readline().strip("{}\n").split(", ")))

### set operations ###
union = set_A.union(set_B)
intersection = set_A.intersection(set_B)
A_diff = set_A.difference(set_B)
B_diff = set_B.difference(set_A)
A_comp = set_U.difference(set_A)
B_comp = set_U.difference(set_B)

### print ###
print(union)
print(intersection)
print(A_diff)
print(B_diff)
print(A_comp)
print(B_comp)
