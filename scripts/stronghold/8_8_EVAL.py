# 8_8_EVAL
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_eval.txt", "r") as f:
    n = int(f.readline().rstrip())
    s = f.readline().rstrip()
    gc_contents = [float(x) for x in f.readline().rstrip().split()]


print(n)
print(s)
print(gc_contents)
