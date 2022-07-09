from math import comb
from urllib.request import urlopen


def read_seq(f):
    # reading file with only a single sequence, any sequence possible
    seq = []
    for line in f:
        seq.append(line.rstrip())
    return ''.join(seq)


def parse_fasta(f):
    # parsing fasta file, returns iteratable object
    name, seq = None, []
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, "".join(seq)


def parse_uniport(url):
    read = urlopen(url).read()  # returns bytes str (b'')
    read_li = str(read).split(r"\n")  # [b'name, ... "'"]
    name = read_li.pop(0)[2:]  # retrieve and cleanse header
    seq = ""
    for s in read_li[:-1]:  # without last element (')
        seq += s
    return name, seq


def mendel_prob(dom: float, rec: float, phen: str) -> float:
    """
    calculate probability of resulting generation's phenotype when dom, rec given
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
    functions recursively
    :param gen: targeted generation
    :param litter: number of offsprings a mature rabbit produces
    :param init: initial number of rabbit pairs
    :param mem: cached list for faster calculation
    :return: calculated number of rabbits for a targeted generation
    """
    if (mem is None) and (gen >= 3):
        mem = [0] + [init] * 2 + [0] * (gen - 2)
    if gen <= 2:
        return mem[gen]
    else:  # gen >= 3
        if mem[gen] == 0:
            mem[gen] = fib_rabbits(gen - 1, litter, init, mem) + litter * fib_rabbits(gen - 2, litter, init, mem)
            return mem[gen]
        else:  # mem[gen] != 0:
            return mem[gen]


def fib_death(gen, lifespan, litter=1, init=1, mem: list = None) -> tuple:
    """
    calculate the number of newborns and mature rabbits following a altered fibonacci method with a lifespan
    functions recursively
    :param gen: targeted generation
    :param lifespan: the lifespan of a rabbit pair
    :param litter: number of offsprings a mature rabbit produces
    :param init: initial number of rabbit pairs
    :param mem: cached list for faster calculation, contains the number of newborns/matures separately
    """
    if mem is None:
        mem = [(0, 0)] + [(init, 0)] + [(0, 0)] * (gen-1)

    if gen == 1:
        return mem[1]
    else:  # gen >= 2
        if sum(mem[gen]) == 0:
            if gen <= lifespan:
                mem[gen] = (fib_death(gen-1, lifespan, litter, init, mem)[1] * litter,
                            sum(fib_death(gen-1, lifespan, litter, init, mem)))
            else:  # gen > lifespan
                mem[gen] = (fib_death(gen - 1, lifespan, litter, init, mem)[1] * litter,
                            sum(fib_death(gen - 1, lifespan, litter, init, mem)) - fib_death(gen-lifespan, lifespan, litter, init, mem)[0])
            return mem[gen]
        else:  # cached
            return mem[gen]


def substring(sup, sub, base: int = 1) -> list:
    """
    finds indices where substrings appear in a sup_string
    function initially checks if first and last nucleotide of the substring appears within the superstring
        with the same distance seen in a substring
    :param sup: super-string to be searched against
    :param sub: sub-string to be searched of
    :param base: 1-based numbering or 0-based numbering
    :return: list containing the start indices where sub appears in sup
    """
    idxs = []
    sup_len, sub_len = map(len, [sup, sub])
    for idx in range(sup_len - sub_len):
        if (sup[idx] == sub[0]) and (sup[idx + sub_len - 1]):
            if sup[idx:idx + sub_len] == sub:
                idxs.append(idx + base)
    return idxs


def binomial_sum(p, n, coef: bool = False) -> float:
    """
    Returns a binomial theorem sum
    :param p: probability p
    :param n: integer n, total
    :param coef: default False, configured when additional multiplication is required for each term
    :return: binomial theorem sum
    """
    q = 1 - p
    res = 0
    if not coef:
        for k in range(0, n + 1):
            res += comb(n, k) * (p ** k) * (q ** (n - k))
    else:  # coef == True
        for k in range(0, n + 1):
            res += k * comb(n, k) * (p ** k) * (q ** (n - k))
    return res


def binomial_sum_limit(p, n, limit, coef: bool = False) -> float:
    """
    binomial_sum with a right end limit
    :param p: probability p
    :param n: integer n, total
    :param limit: right end limit (from 0 to limit)
    :param coef: default false, configured when additional multiplication is required for each term
    :return: binomial theorem sum
    """
    q = 1 - p
    res = 0
    if not coef:
        for k in range(0, limit + 1):
            res += comb(n, k) * (p ** k) * (q ** (n - k))
    else:  # coef == True
        for k in range(0, limit + 1):
            res += k * comb(n, k) * (p ** k) * (q ** (n - k))
    return res


def overlap(seq1, seq2, n: int) -> bool:
    """
    return whether the length n suffix of seq1 matches length n prefix of seq2
    :param seq1: sequence which its suffix is matched against
    :param seq2: sequence which its prefix is matched against
    :param n: overlap length
    """
    return seq1[len(seq1)-n:] == seq2[:n]


# Nucleotides
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

    def gc_content(self, percent=True) -> float:
        gc = 0
        for nt in self.seq:
            if (nt == "G") or (nt == "C"):
                gc += 1
        if percent:
            return gc / len(self.seq) * 100
        else:  # percent = False
            return gc / len(self.seq)

    def translate(self) -> str:
        codon_table = {"UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
                       "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
                       "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
                       "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
                       "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
                       "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
                       "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
                       "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
                       "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
                       "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
                       "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
                       "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E",
                       "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
                       "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
                       "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
                       "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}
        aa_seq = ''
        for idx in range(0, len(self.seq), 3):
            codon = self.seq[idx:idx + 3]
            aa = codon_table[codon]
            if aa == "Stop":
                break
            else:  # not an stop codon
                aa_seq += aa
        return aa_seq


# Proteins
class proteins:

    def __init__(self, seq):
        self.seq = seq

    def reverse_translate_var(self, end=True) -> int:
        reverse_table = {"A": ("GCU", "GCC", "GCA", "GCG"), "C": ("UGU", "UGC"), "D": ("GAU", "GAC"),
                         "E": ("GAA", "GAG"), "F": ("UUU", "UUC"), "G": ("GGU", "GGC", "GGA", "GGG"),
                         "H": ("CAU", "CAC"), "I": ("AUU", "AUC", "AUA"), "K": ("AAA", "AAG"),
                         "L": ("UUA", "UUG", "CUU", "CUC", "CUA", "CUG"), "M": ("AUG",), "N": ("AAU", "AAC"),
                         "P": ("CCU", "CCC", "CCA", "CCG"), "Q": ("CAA", "CAG"),
                         "R": ("CGU", "CGC", "CGU", "CGC", "CGA", "CGG"),
                         "S": ("UCU", "UCC", "UCA", "UCG", "AGU", "AGC"), "T": ("ACU", "ACC", "ACA", "ACG"),
                         "V": ("GUU", "GUC", "GUA", "GUG"), "W": ("UGG",), "Y": ("UAU", "UAC"),
                         "stop": ("UAA", "UAG", "UGA")}
        var = 1
        if end:  # consider stop codon variability
            var *= 3
        for aa in self.seq:
            var *= len(reverse_table[aa])
        return var

    def weight(self) -> float:
        # monoisotopic mass table
        mass_table = {"A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259, "F": 147.06841, "G": 57.02146,
                      "H": 137.05891, "I": 113.08406, "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
                      "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203, "T": 101.04768, "V": 99.06841,
                      "W": 186.07931, "Y": 163.06333}
        mass = 0
        for aa in self.seq:
            mass += mass_table[aa]
        return mass
