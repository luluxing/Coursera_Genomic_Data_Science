from Bio.Seq import Seq

def read_fasta(file):
	from Bio import SeqIO
	handle = open(file, "rU") 
	# 'rU' means open for reading using universal readline mode 
	seqs = {}
	for record in SeqIO.parse(handle, "fasta"):
		name, sequence = record.id, str(record.seq)
		if name not in seqs:
			seqs[name] = sequence
		else:
			print "Error! Same identifications in one file"
	
	handle.close()
	return seqs

def translate(s, frame=1):
	# if frame != 2 and frame != 3:
	# 	return "Error! Illegal frame number!"
	stops = ['TAA', 'TAG', 'TGA']
	newseq = s[frame-1:]
	i = 0
	neworf = False
	maxorflen = 0
	tmp = 0
	while i < len(newseq):
		if neworf and i + 3 <= len(newseq):
			tmp += 3
			if newseq[i:i+3] in stops:
				if tmp >= maxorflen:
					maxorflen = tmp
					maxini = ini
				tmp = 0
				neworf = False
		elif not neworf and i + 3 <= len(newseq):
			if newseq[i:i+3] == 'ATG':
				ini = i
				tmp += 3
				neworf = True
		i += 3
		
	return maxorflen

def construct_repeat_dic():
	repdic = {}
	from itertools import product    
	c = [(''.join(codon)) for codon in product('ATCG', repeat=6)]
	for k in c:
		repdic[k] = 0
	return repdic

def find_repeats(s,replen,allrep):
	# allrep = construct_repeat_dic()
	i = 0
	while i < len(s):
		if i + replen <= len(s): 
			tmp =  s[i:i+replen]
			if tmp in allrep:
				allrep[tmp] += 1
			else:
				allrep[tmp] = 1
		i += 1
	return allrep
	# return allrep[max(allrep, key=lambda i: allrep[i])]


seqs = read_fasta("dna2.fasta")

# # Question 1
# print len(seqs)

# # Question 2&3
# maxseqlen = 0
# minseqlen = 100000
# for k in seqs:
# 	if len(seqs[k]) >= maxseqlen:
# 		maxseqlen = len(seqs[k])
# 	if len(seqs[k]) <= minseqlen:
# 		minseqlen = len(seqs[k])
# print maxseqlen, minseqlen

# Question 4
# for k in seqs:
# 	print k
# 	print translate(seqs[k], 2)

# Question 5
# for k in seqs:
# 	a = translate(seqs[k], 1)
# 	b = translate(seqs[k], 2)
# 	c = translate(seqs[k], 3)
# 	print k,max(a,b,c)

# Question 8
allrep = {}
import operator
for k in seqs:
	i = 0
	s = seqs[k]
	replen = 7
	while i < len(s):
		if i + replen <= len(s): 
			tmp =  s[i:i+replen]
			if tmp in allrep:
				allrep[tmp] += 1
			else:
				allrep[tmp] = 1
		i += 1

kys = ['TGCGCGC','CGCGCCG','CATCGCC','GCGGCCG']
for k in kys:
	print allrep[k]

# sorted_allrep = sorted(allrep.items(), key=operator.itemgetter(1))
# print sorted_allrep[44711:],len(sorted_allrep)
# print allrep[max(allrep, key=lambda i: allrep[i])]