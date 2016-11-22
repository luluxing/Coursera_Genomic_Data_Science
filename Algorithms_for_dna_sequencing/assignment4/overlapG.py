from readFq import readFastq
from overlap import overlap

def kmersetDic(reads, k):
	setdic = {}
	for read in reads:
		for i in range(0, len(read)-k+1):
			if read[i:i+k] not in setdic:
				setdic[read[i:i+k]] = set()
			setdic[read[i:i+k]].add(read)
	return setdic

def overlap_all_pairs(reads, k):
	kmerset = kmersetDic(reads, k)

	overlap_g = {}
	
	for read in reads:
		read_sfx = read[len(read)-k:]
		for readb in kmerset[read_sfx]:
			if readb == read:
				continue
			if overlap(read, readb, k) != 0:
				if read not in overlap_g:
					overlap_g[read] = []
				overlap_g[read] += [readb]

	return overlap_g




