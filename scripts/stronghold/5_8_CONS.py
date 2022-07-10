# 5_8_CONS
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_cons.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))

### create profile matrix ###
profile = rs.profile_matrix(fastd, consensus=False, nt_order="ACGT")
consensus = rs.profile_matrix(fastd, consensus=True, nt_order="ACGT")

### print results ###
print(consensus)
for nt in "ACGT":
    print(" ".join([f"{nt}:"] + [str(x) for x in profile.loc[nt, :]]))
