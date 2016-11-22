from naive_with_counts import naive_with_counts
from bm_with_counts import boyer_moore_with_counts
from bm_preproc import BoyerMoore
from kmer_index import Index
from subseq_index import SubseqIndex

def readFasta(filename):
	with open(filename, 'r') as f:
		sequence = ''
		for line in f:
			if line[0] == '>':
				continue
			sequence += line.rstrip()
	f.close()
	return sequence

def query_2mm(p, t, hit, pg_num):
	if pg_num == 1:
		if hit + 16 >= len(t):
			return False
		p_new = p[8:]
		t_new = t[hit+8:hit+16]
	elif pg_num == 2:
		if hit - 8 < 0 or hit + 16 >= len(t):
			return False
		p_new = p[:8] + p[16:]
		t_new = t[hit-8:hit] + t[hit+8:hit+16]
	else:
		if hit - 16 < 0:
			return False
		p_new = p[:16]
		t_new = t[hit-16:hit]

	if len(p_new) != len(t_new):
		return False
	mismatch = 0
	for i in range(len(p_new)):
		if p_new[i] != t_new[i]:
			mismatch += 1
		if mismatch > 2:
			return False
	return True

def kmer_2mm(p, t):
	hits_num = 0
	if len(p) != 24:
		return 'ERROR! Length not 24nt!'
	t_8mer = Index(t, 8)
	hit_2mm = []
	for i in range(1, 4):
		hits = t_8mer.query(p[(i-1)*8:i*8])
		for hit in hits:
			hits_num += 1
			if query_2mm(p, t, hit, i) and (hit - (i-1)*8) not in hit_2mm:
				hit_2mm.append((hit - (i-1)*8))
	return hit_2mm, hits_num


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatch = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatch += 1
            if mismatch > 2:
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def query_subseq(p, t, subseq_ind):
	hits_1 = subseq_ind.query(p[0:])
	hits_2 = subseq_ind.query(p[1:])
	hits_3 = subseq_ind.query(p[2:])
	hit_num = len(hits_3+hits_2+hits_1)
	return hit_num


hg38_chr1 = readFasta('chr1.GRCh38.excerpt.fasta')

#q1
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, hg38_chr1)
print num_alignments # 799954

#q2
print num_character_comparisons # 984143

# #q3
# p_bm = BoyerMoore(p)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, hg38_chr1)
# print num_alignments # 127974

#q4
# p = 'GGCGCGGTGGCTCACGCCTGTAAT'
# match, hits_num =  kmer_2mm(p, hg38_chr1)
# match, hits_num = naive_2mm(p, hg38_chr1)
# print len(match) # 19

#q5
# print hits_num # 90

#q6
# subseq_ind = SubseqIndex(hg38_chr1, 8, 3)
# print query_subseq(p, hg38_chr1, subseq_ind) # 79
