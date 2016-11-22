def naive_with_counts(p, t):
    occurrences = []
    aligh_num, chr_comp_num = 0, 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        aligh_num += 1
        for j in range(len(p)):  # loop over characters
            chr_comp_num += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, aligh_num, chr_comp_num

# p = 'word'
# t = 'there would have been a time for such a word'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# # ([40], 41, 46)

# p = 'needle'
# t = 'needle need noodle needle'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# ([0, 19], 20, 35)
