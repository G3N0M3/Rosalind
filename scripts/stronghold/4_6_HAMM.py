# 4_6_HAMM

### read data ###
with open("../../inputs/stronghold/rosalind_hamm.txt", "r") as f:
    s, t = map(lambda x: x.rstrip(), list(f))

### calculate Hamming distance
hamm = 0
for idx in range(len(s)):
    if s[idx] != t[idx]:
        hamm += 1

print(hamm)
