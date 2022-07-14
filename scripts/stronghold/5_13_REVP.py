# 5_13_REVP
import scripts.Rosalind as rs

### read data ###
with open("../../inputs/stronghold/rosalind_revp.txt", "r") as f:
    fastd = rs.make_fastd(rs.parse_fasta(f))
    seq = list(fastd.values())[0]

### check reverse palindrome information ###
# seq = "ATGCGCAT"
len_limit = (4, 12)  # sequence length boundaries for palindromes from~to
pal_info = []
for start in range(0, len(seq) - 1):
    for end in range(start + 1, len(seq)):
        end += 1
        ln = end - start
        if (ln < len_limit[0]) or (ln > len_limit[1]):
            continue  # continue loop to next step if palindrome length out of limit
        # check if palindrome and save info
        pal_seq = seq[start:end]
        if rs.check_reverse_palindrome(pal_seq):
            pal_info.append((start + 1, ln))

### print palindrome information ###
for item in pal_info:
    print(item[0], item[1])
