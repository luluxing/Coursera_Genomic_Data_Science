from readFq import readFastq

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

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

	overlap_g = []
	count = 0
	
	for read in reads:
		read_sfx = read[len(read)-k:]
		flag = False
		for readb in kmerset[read_sfx]:
			if readb == read:
				continue
			if overlap(read, readb, k) != 0:
				overlap_g.append((read, readb))
				if not flag:
					count += 1
					flag = True

	return overlap_g, count



# sequences = readFastq('ERR266411_1.for_asm.fastq')
# print sequences

# reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
# readsoverlap, count = overlap_all_pairs(sequences, 30)
# print len(readsoverlap), count # 904746 7161





