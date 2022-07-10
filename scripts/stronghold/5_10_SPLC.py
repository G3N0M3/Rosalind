# 5_10_SPLC
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_splc.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seqs = list(fastd.values())
    sup = seqs[0]
    introns = seqs[1:]

exons = rs.erase_intron(sup, introns)
exons = rs.nucleotides(exons).transcribe()
aa_seq = rs.nucleotides(exons).translate()
print(aa_seq)
