from math import comb
from urllib.request import urlopen
import pandas as pd

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
               "UAA": "stop", "CAA": "Q", "AAA": "K", "GAA": "E",
               "UAG": "stop", "CAG": "Q", "AAG": "K", "GAG": "E",
               "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
               "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
               "UGA": "stop", "CGA": "R", "AGA": "R", "GGA": "G",
               "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}

reverse_codon_table = {"A": ("GCU", "GCC", "GCA", "GCG"), "C": ("UGU", "UGC"), "D": ("GAU", "GAC"), "E": ("GAA", "GAG"),
                       "F": ("UUU", "UUC"), "G": ("GGU", "GGC", "GGA", "GGG"), "H": ("CAU", "CAC"),
                       "I": ("AUU", "AUC", "AUA"), "K": ("AAA", "AAG"), "L": ("UUA", "UUG", "CUU", "CUC", "CUA", "CUG"),
                       "M": ("AUG",), "N": ("AAU", "AAC"), "P": ("CCU", "CCC", "CCA", "CCG"), "Q": ("CAA", "CAG"),
                       "R": ("CGU", "CGC", "CGU", "CGC", "CGA", "CGG"), "S": ("UCU", "UCC", "UCA", "UCG", "AGU", "AGC"),
                       "T": ("ACU", "ACC", "ACA", "ACG"), "V": ("GUU", "GUC", "GUA", "GUG"), "W": ("UGG",),
                       "Y": ("UAU", "UAC"), "stop": ("UAA", "UAG", "UGA")}

comp_nt = {"A": "T", "T": "A", "G": "C", "C": "G"}

# monoisotopic mass table
mass_table = {"A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259, "F": 147.06841, "G": 57.02146,
              "H": 137.05891, "I": 113.08406, "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
              "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203, "T": 101.04768, "V": 99.06841,
              "W": 186.07931, "Y": 163.06333}


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


def make_fastd(fasta) -> dict:
    fastd = {}
    for name, seq in fasta:
        fastd[name] = seq
    return fastd


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
        mem = [(0, 0)] + [(init, 0)] + [(0, 0)] * (gen - 1)

    if gen == 1:
        return mem[1]
    else:  # gen >= 2
        if sum(mem[gen]) == 0:
            if gen <= lifespan:
                mem[gen] = (fib_death(gen - 1, lifespan, litter, init, mem)[1] * litter,
                            sum(fib_death(gen - 1, lifespan, litter, init, mem)))
            else:  # gen > lifespan
                mem[gen] = (fib_death(gen - 1, lifespan, litter, init, mem)[1] * litter,
                            sum(fib_death(gen - 1, lifespan, litter, init, mem)) -
                            fib_death(gen - lifespan, lifespan, litter, init, mem)[0])
            return mem[gen]
        else:  # cached
            return mem[gen]


def substring(sup, sub) -> bool:
    """
    return if sub is a substring of sup
    :param sup: sequence to be searched upon
    :param sub: sequence to be searched of
    :return: boolean
    """
    match = False
    sup_len, sub_len = map(len, [sup, sub])
    for idx in range(sup_len - sub_len):
        if (sup[idx] == sub[0]) == (sup[idx + sub_len - 1] == sub[-1]):
            if sup[idx:idx + sub_len] == sub:
                match = True
    return match


def substring_idx(sup, sub, base: int = 1) -> list:
    """
    finds indices where substrings appear in a sup_string
    function initially checks if first and last nucleotide of the substring appears within the superstring
        with the same distance seen in a substring
    :param sup: super-string to be searched upon
    :param sub: sub-string to be searched of
    :param base: 1-based numbering or 0-based numbering
    :return: list containing the start indices where sub appears in sup
    """
    idxs = []
    sup_len, sub_len = map(len, [sup, sub])
    for idx in range(sup_len - sub_len):
        if (sup[idx] == sub[0]) and (sup[idx + sub_len - 1] == sub[-1]):
            if sup[idx:idx + sub_len] == sub:
                idxs.append(idx + base)
    return idxs


def shared_motif(seqs: list) -> str:
    """
    returns the longest shared motif of a given list of sequences
    :param seqs: list of sequences
    :return: longest shared motif
    """
    # use the shortest sequence as the sequence to get the shared motif
    seqs = sorted(seqs, key=len)
    sub_seq = seqs.pop(0)
    # finding longest shared motif
    common = ""
    for ln in range(len(sub_seq), 1, -1):
        for start_idx in range(0, len(sub_seq) - ln + 1):
            sub = sub_seq[start_idx:start_idx + ln]
            seq, match = "", False
            for seq in seqs:
                match = substring(seq, sub)
                if not match:  # if substring does not exits
                    break  # break to assign new seq
            if match and (seq == seqs[-1]):
                common = sub
                break
        if common != "":  # if common substring is assigned
            break
    return common


def erase_intron(sup, introns: list) -> str:
    """
    eliminates intron part from given sup-string
    function erase sequences given in "introns" from 0 to end (order-specific workflow)
    :param sup: sup-string containing both exons and introns
    :param introns: list of sequences to be erased from sup-string
    :return: intron-erased string
    """
    for intron in introns:
        # get indexes
        idxs = substring_idx(sup, intron, base=0)
        # erase intron if exists
        if len(idxs) != 0:
            idx = idxs[0]
            # if intron is at end of sup
            if idx + len(intron) == len(sup):
                sup = sup[:idx + 1]
            else:  # intron in the middle of sup
                sup = sup[:idx] + sup[idx + len(intron):]

    return sup


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
    return seq1[len(seq1) - n:] == seq2[:n]


def fragment_assembly(seq_list) -> str:
    """
    returns shortest superstring from a list of sequence with equal length
    :param seq_list: list of sequences with equal length
    :return: shortest superstring
    """
    res = seq_list.pop(0)
    n = len(res)
    while len(seq_list) > 0:
        for seq_idx in range(len(seq_list)):
            seq = seq_list[seq_idx]
            match = False
            for ln in range((n - 1), int(n / 2) - 1, -1):
                # "res + seq" match
                if overlap(res, seq, n=ln):
                    match = True
                    res = res + seq[ln:]
                    del seq_list[seq_idx]
                    break
                # "seq + res" match
                if overlap(seq, res, n=ln):
                    match = True
                    res = seq + res[ln:]
                    del seq_list[seq_idx]
                    break
            # if there was a match, continue whole progress with new res
            if match:
                break  # break -> still in while loop
    return res


def profile_matrix(fastd, consensus: bool = False, nt_order: str = "ACGT"):
    """
    calculate profile matrix and consensus string from a given fastd
    :param fastd: fastd containing names as keys and sequences as values
    :param consensus: if function should return profile matrix or consensus string
    :param nt_order: order of nucleotides in profile matrix
    :return: profile matrix in pd.DataFrame or consensus string, according to the consensus parameter
    """
    seqs = list(fastd.values())
    string_len = len(seqs[0])
    profile = pd.DataFrame(0, index=list(nt_order),
                           columns=[x + 1 for x in range(string_len)])
    for seq in seqs:
        for idx in range(string_len):
            nt = seq[idx]
            profile.loc[nt, idx + 1] += 1
    consensus_str = "".join(profile.idxmax())

    # return
    if consensus:
        return consensus_str
    else:  # return profile matrix
        return profile


def check_reverse_palindrome(seq) -> bool:
    # configure iteration range
    n = len(seq)
    if n % 2 == 0:  # length of sequence is even
        rng = int(n / 2)
    else:  # length of sequence is odd -> can not be a reverse palindrome
        pal = False
        return pal

    # check if reverse palindrome
    for idx in range(rng):
        pal = (comp_nt[seq[idx]] == seq[n - idx - 1])
        if not pal:  # is there a mismatch in corresponding nucleotides?
            break

    return pal


def permutate_ascend(to_list: list, from_list: list) -> list:
    """
    function made to be used for function "spliced_motif"
    appends numbers from from_list to to_list in a ascending pattern
    :param to_list: list to be appended to
    :param from_list: list to choose number to be appended from
    :return: list of lists with numbers in a ascending pattern
    """
    res = []
    for to_item in to_list:
        for from_num in from_list:
            if to_item[-1] < from_num:
                res.append(to_item + [from_num])
    return res


def search_spliced_motif(sup, motif, case, base: int = 1) -> list:
    """
    returns list of spliced motifs
    !! theoretically works, but takes too much time, even after refining sub_array
    -> solved by early cutting in permutation loop
    :param sup: sup-string
    :param motif: motif to be search in sup-string
    :param case: number of cases to be returned
    # -1 to retrieve all spliced motifs, but unrecommended for long motifs or long sequences
    :param base: 1-based numbering or 0-based numbering
    # all calculations are done by 0-based numbering
    # if 1-base numbering is to be returned, indexes are modified right before actual return
    :return: list of spliced motifs
    """
    ## parse indexes in sup
    a_idx, t_idx, g_idx, c_idx = [], [], [], []
    for nt_idx in range(len(sup)):
        nt = sup[nt_idx]
        if nt in ("A", "T"):
            if nt == "A":
                a_idx.append(nt_idx)
            else:  # nt == "T"
                t_idx.append(nt_idx)
        else:  # nt in ("G", "C")
            if nt == "G":
                g_idx.append(nt_idx)
            else:  # nt == "C"
                c_idx.append(nt_idx)
    nt_idx = {"A": a_idx, "T": t_idx, "G": g_idx, "C": c_idx}

    ## make motif array
    motif_array = []
    for nt in motif:
        motif_array.append(nt_idx[nt][:])  # [:] for shallow copy?

    ## refine sub_array
    # refine elements in the front of each row
    for row in range(len(motif_array) - 1):
        erase_crit = motif_array[row][0]
        while motif_array[row + 1][0] <= erase_crit:
            del motif_array[row + 1][0]
    # refine elements in the end of each row
    for row in range(len(motif_array) - 1, 0, -1):
        erase_crit = motif_array[row][-1]
        while motif_array[row - 1][-1] >= erase_crit:
            del motif_array[row - 1][-1]

    ## case calculation for total retrieve
    if case == -1:
        case = 1
        for mtf in motif_array:
            case *= len(motif)

    ## permutate to get every spliced motifs
    to_list = [[x] for x in motif_array.pop(0)]
    for from_list in motif_array:
        to_list = permutate_ascend(to_list=to_list, from_list=from_list)
        # the following [:case] controls number of spliced motifs to be returned
        to_list = to_list[:case]  # critical for time-consumption control

    ## final to_list is result to return
    res = to_list

    ## apply base numbering
    if base:  # if 1-base numbering
        for idx in range(len(res)):
            res[idx] = [sum(x) for x in zip(res[idx], [1] * len(res[idx]))]

    return res


def point_mutation(nt1, nt2) -> str:
    """
    return "same" (no mutation), "transition", "transversion"
    :param nt1: nucleotide from original sequence
    :param nt2: nucleotide from mutated sequence
    """
    nt_type = {"A": "purine", "T": "pyrimidine", "G": "purine", "C": "pyrimidine"}
    if nt1 == nt2:
        return "same"
    else:  # point mutation
        if nt_type[nt1] == nt_type[nt2]:
            return "transition"
        else:
            return "transversion"


def longest_pattern_subsequence(seq_list: list, ascending: bool = True) -> str:
    len_list = []
    n = len(seq_list)

    for num in seq_list:
        for idx in range(n):  # from sequence length of 1 to n
            if len(len_list) < (idx + 1):
                len_list.append([])
            _seq = len_list[idx]
            if idx == 0:
                if not _seq:  # empty
                    len_list[idx] = [num]
                    break  # to next number
                else:  # not empty
                    if ascending:
                        if _seq[-1] > num:
                            len_list[idx] = [num]
                            break
                        else:
                            continue  # to next length
                    else:  # descending
                        if _seq[-1] < num:
                            len_list[idx] = [num]
                            break
                        else:
                            continue  # to next length
            else:  # length bigger than 1
                if not _seq:  # empty
                    if ascending:
                        if len_list[idx - 1][-1] < num:
                            len_list[idx] = len_list[idx - 1] + [num]
                            break  # to next number
                        else:
                            pass  # CHECK!!!
                    else:  # descending
                        if len_list[idx - 1][-1] > num:
                            len_list[idx] = len_list[idx - 1] + [num]
                            break
                        else:
                            pass  # CHECK!!!
                else:  # not empty
                    if ascending:
                        if _seq[-1] > num:
                            len_list[idx] = len_list[idx - 1] + [num]
                            break  # to next number
                        else:
                            continue  # to next length
                    else:  # descending
                        if _seq[-1] < num:
                            len_list[idx] = len_list[idx - 1] + [num]
                            break
                        else:
                            continue

    return " ".join(map(str, len_list[-1]))


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
        res = ''
        for nt in self.seq:
            res += comp_nt[nt]
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
        aa_seq = ''
        for idx in range(0, len(self.seq), 3):
            codon = self.seq[idx:idx + 3]
            aa = codon_table[codon]
            if aa == "stop":
                break
            else:  # not an stop codon
                aa_seq += aa
        return aa_seq

    def orf(self, rev_comp: bool = True) -> list:
        """
        return list with all possible
        :param rev_comp:
        :return:
        """
        # define targeted sequences
        seq = nucleotides(self.seq)
        seqs = [seq.transcribe()]
        if rev_comp:
            seq_comp = nucleotides(seq.complement(reverse=True)).transcribe()
            seqs.append(seq_comp)

        # get ORFs
        orf_li = []
        for sq in seqs:
            for base in [0, 1, 2]:
                for start_idx in range(base, len(sq) - 3 + 1, 3):
                    aa_seq = ""
                    codon = sq[start_idx:start_idx + 3]
                    aa = codon_table[codon]
                    if (aa_seq == "") and (aa != "M"):  # start codon not Met
                        continue
                    # code below proceeds when start codon is Met
                    aa_seq += aa
                    for idx in range(start_idx + 3, len(sq) - 3, 3):
                        codon = sq[idx:idx + 3]
                        aa = codon_table[codon]
                        # end sequence of append amino acid to sequence
                        if aa == "stop":
                            orf_li.append(aa_seq)
                            break
                        else:
                            # last possible codon, but not a stop codon
                            # -> erase sequence and proceed to next start index
                            if idx + 3 > len(sq):
                                break
                            # not a stop codon, not the last possible codon
                            aa_seq += aa
                            continue

        orf_li = list(set(orf_li))

        return orf_li

    def check_reverse_palindrome(self):
        return check_reverse_palindrome(self.seq)


# Proteins
class proteins:

    def __init__(self, seq):
        self.seq = seq

    def reverse_translate_var(self, end=True) -> int:
        var = 1
        if end:  # consider stop codon variability
            var *= 3
        for aa in self.seq:
            var *= len(reverse_codon_table[aa])
        return var

    def weight(self) -> float:
        mass = 0
        for aa in self.seq:
            mass += mass_table[aa]
        return mass

    def n_gly_motif(self, base: int = 1) -> list:
        """
        return indexes for N-glycosylation motif (default 1-based)
        motif -> N{P}[ST]{P}
        [XY] -> X or Y / {X} -> anything but X
        :param base:defines if 1-based or 0-based for index numbering
        :return list containing indexes for where the motif starts
        """
        res = []
        for idx in range(len(self.seq) - 4):
            if self.seq[idx] == "N":
                if self.seq[idx + 1] != "P":
                    if (self.seq[idx + 2] == "S") or (self.seq[idx + 2] == "T"):
                        if self.seq[idx + 3] != "P":
                            res.append(idx + base)
        return res
