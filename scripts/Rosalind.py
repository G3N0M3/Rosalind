# Nucleotides

def read_seq(f):
    seq = []
    for line in f:
        seq.append(line.rstrip())
    return ''.join(seq)


def parse_fasta(f):
    name, seq = None, []
    for line in f:
        line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
        else:
            seq.append(line)
    if name:
        yield name, "".join(seq)


class nucleotides:

    def __init__(self, seq):
        self.seq = seq

    def count(self) -> dict:
        # Base count
        counts = {"A": 0, "T": 0, "G": 0, "C": 0, "U": 0}
        for nt in self.seq:
            counts[nt] += 1

        return counts

    def transcribe(self) -> str:
        return self.seq.replace("T", "U")

    def complement(self, reverse=False) -> str:
        comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
        res = ''
        for nt in self.seq:
            res += comp[nt]

        if reverse:
            res = res[::-1]

        return res
