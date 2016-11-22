from bm_preproc import BoyerMoore

def boyer_moore_with_counts(p, p_bm, t):
	i = 0
	occurrences = []
	align_num, chr_comp_num = 0, 0
	while i < len(t) - len(p) + 1:
		align_num += 1
		shift = 1
		mismatched = False
		for j in range(len(p)-1, -1, -1):
			chr_comp_num += 1
			if p[j] != t[i+j]:
				skip_bc = p_bm.bad_character_rule(j, t[i+j])
				skip_gs = p_bm.good_suffix_rule(j)
				shift = max(shift, skip_bc, skip_gs)
				mismatched = True
				break
		if not mismatched:
			occurrences.append(i)
			skip_gs = p_bm.match_skip()
			shift = max(shift, skip_gs)
		i += shift
	return occurrences, align_num, chr_comp_num

# p = 'word'
# t = 'there would have been a time for such a word'
# lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)
# ([40], 12, 15)

# p = 'needle'
# t = 'needle need noodle needle'
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)
# ([0, 19], 5, 18)

