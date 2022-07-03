# 4_3_GC
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_gc.txt", "r") as f:
    fastl = list(rs.parse_fasta(f))

### calculate GC content
gc_max = -1
name_max = ''
for name, seq in fastl:
    gc = rs.nucleotides(seq).gc_content(percent=True)
    if gc > gc_max:
        gc_max = gc
        name_max = name[1:]

print(name_max)
print(gc_max)
