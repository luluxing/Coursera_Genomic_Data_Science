from readFq import readFastq
import overlapG
from overlap import overlap
import itertools, random

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def pick_max_overlap(reads, kmerset, k):
    # bestreads = []
    bestreads = ()
    max_overlap = 0
    for read in reads:
        read_sfx = read[len(read)-k:]
        for readb in kmerset[read_sfx]:
            if readb == read:
                continue
            olen = overlap(read, readb, k)
            if olen > max_overlap:
                # bestreads = [(read, readb)]
                bestreads = (read, readb)
                max_overlap = olen
    return bestreads, max_overlap
    #         elif olen == max_overlap:
    #             bestreads.append((read, readb))
    # if bestreads != []:
    #     return random.choice(bestreads), max_overlap
    # return [], 0

def greedy_scs(reads, k):
    kmerset = overlapG.kmersetDic(reads, k)
    pair, maxlen = pick_max_overlap(reads, kmerset, k)
    while maxlen > 0:
        reada, readb = pair[0], pair[1]
        new_read = reada + readb[maxlen:]

        reads.remove(reada)
        reads.remove(readb)
        reads.append(new_read)

        for i in range(0, len(new_read)-k+1):
            kmerset[new_read[i:i+k]].discard(reada)
            kmerset[new_read[i:i+k]].discard(readb)
            kmerset[new_read[i:i+k]].add(new_read)
        
        pair, maxlen = pick_max_overlap(reads, kmerset, k)

    return ''.join(reads)



virusfq = readFastq('ads1_week4_reads.fq')
assembled = greedy_scs(virusfq, 40)
# assembled = dfs(debruijn(virusfq, 50))
print len(assembled), assembled.count('A'), assembled.count('T')











