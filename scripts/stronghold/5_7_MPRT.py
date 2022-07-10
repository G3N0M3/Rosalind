# 5_7_MPRT
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_mprt.txt", "r") as f:
    access_ids = [x.rstrip() for x in f]

### parse uniport data
fastd = {}
for prt_id in access_ids:
    prt_url = ''.join(["http://www.uniprot.org/uniprot/",
                       prt_id,
                       ".fasta"])
    name, seq = rs.parse_uniport(prt_url)
    name = name.split("|")[1]  # retrieve protein name from header
    fastd[name] = seq

### find motif indexes ###
for prt_id in access_ids:
    prt_seq = rs.proteins(fastd[prt_id])
    motif_idxs = prt_seq.n_gly_motif(base=1)
    if len(motif_idxs) != 0:
        print(prt_id)
        print(" ".join([str(x) for x in motif_idxs]))
