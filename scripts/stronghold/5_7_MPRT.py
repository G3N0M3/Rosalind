# 5_7_MPRT
import scripts.Rosalind as rs
from urllib.request import urlopen

### read data ###
with open("../../inputs/stronghold/rosalind_mprt.txt", "r") as f:
    access_ids = [x.rstrip() for x in f]

###
fastd = {}
for prt_id in access_ids:
    prt_url = ''.join(["http://www.uniprot.org/uniprot/",
                       prt_id,
                       ".fasta"])
    name, seq = rs.parse_uniport(prt_url)
    name = name.split("|")[1]  # retrieve protein name from header
    fastd[name] = seq

print(fastd)
