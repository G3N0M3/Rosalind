# Nucleotides

def read_seq(f):
    # reading file with only a single sequence
    seq = []
    for line in f:
        seq.append(line.rstrip())
    return ''.join(seq)


def parse_fasta(f):
    # parsing fasta file, returns iteratable object
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


def mendel_prob(dom: float, rec: float, phen: str) -> float:
    """
    calculate probability of resulting generation's phenotype
    :param dom: ratio of dominant allele
    :param rec: ratio of recessive allele
    :param phen: targeted phenotype of offsprings to be calculated
    :return: probability of targeted phenotype shown in offspring generation
    """
    if phen == "dominant":
        prob = dom ** 2 + 2 * dom * rec
    elif phen == "recessive":
        prob = rec ** 2
    else:
        raise ValueError("Choose btw 'dominant' or 'recessive'")
    return prob


def fib_rabbits(gen, litter=1, init=1, mem: list = None):
    """
    calculate the number of rabbits following a altered fibonacci method
    :param gen: targeted generation which
    :param litter: number of offsprings a mature rabbit produces
    :param init: initial number of rabbits
    :param mem: cached list for faster calculation
    :return: calculated number of rabbits for a targeted generation
    """
    if (mem is None) and (gen >= 3):
        mem = [0] + [litter]*2 + [0]*(gen-2)
    if gen <= 2:
        return init
    else:  # gen >= 3
        if mem[gen] == 0:
            mem[gen] = fib_rabbits(gen - 1, litter, init, mem) + litter * fib_rabbits(gen - 2, litter, init, mem)
            return mem[gen]
        else:  # mem[gen] != 0:
            return mem[gen]


class nucleotides:
    # class used for handling nucleotide sequences

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
        # reverse complement
        if reverse:
            res = res[::-1]
        return res
