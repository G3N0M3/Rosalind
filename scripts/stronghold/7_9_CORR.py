# 7_9_KMER
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_corr.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seqs = list(fastd.values())

### count appearance of sequences ###
seq_counts = {}
for seq in seqs:
    if seq in seq_counts:
        seq_counts[seq] += 1
    else:  # seq is new
        # check if reverse-complement exists
        rev_comp = rs.reverse_complement(seq)
        if rev_comp in seq_counts:
            seq_counts[rev_comp] += 1
        else:  # if seq is "completely" new
            seq_counts[seq] = 1

### identify correct/incorrect reads ###
"""
correct reads appear at least twice, 
whereas incorrect reads appear only once
"""
corrects, incorrects = [], []
for seq in seq_counts:
    if seq_counts[seq] == 1:  # read is incorrect
        incorrects.append(seq)
    else:
        corrects.append(seq)

### calculate hamming distance and match incorrect reads with correct reads ###
modification_match = []
for incorrect_read in incorrects:
    for correct_read in corrects:
        # find hamming distance == 1 match
        if rs.hamm_distance(incorrect_read, correct_read) == 1:
            modification_match.append(f"{incorrect_read}->{correct_read}")
            break
    # if no match found
    revComp_incorrect = rs.reverse_complement(incorrect_read)
    for correct_read in corrects:
        if rs.hamm_distance(revComp_incorrect, correct_read) == 1:
            revComp_correct = rs.reverse_complement(correct_read)
            modification_match.append(f"{incorrect_read}->{revComp_correct}")  ###
            break

### print ###
for item in modification_match:
    print(item)
