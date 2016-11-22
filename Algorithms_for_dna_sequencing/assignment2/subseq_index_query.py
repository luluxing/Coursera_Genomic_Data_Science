from subseq_index import SubseqIndex


def query_subseq(p, t, subseq_ind):
	hits_1 = subseq_ind.query(p[0:])
	hits_2 = subseq_ind.query(p[1:])
	hits_3 = subseq_ind.query(p[2:])
	hit_num = len(hits_3+hits_2+hits_1)
	print hits_1,hits_2,hits_3

t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = query_subseq(p, t, subseq_ind)
print(occurrences)
# [0, 14]
print(num_index_hits)
# 6