# 6_7_TRAN
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_tran.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))

seqs = list(fastd.values())
original = seqs[0]
mutated = seqs[1]

transition = 0
transversion = 0
for idx in range(len(original)):
    ori_nt = original[idx]
    mut_nt = mutated[idx]

    point_status = rs.point_mutation(ori_nt, mut_nt)
    if point_status == "transition":
        transition += 1
    elif point_status == "transversion":
        transversion += 1

ratio = transition / transversion

print(ratio)
