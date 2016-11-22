def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities



def naive_with_rc(p, t):
    occurrences = naive(p, t)
    occurrences_rc = naive(reverseComplement(p), t)
    for offset in occurrences_rc:
        if offset not in occurrences:
            occurrences.append(offset)
    return occurrences

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

phix_genome = readGenome('lambda_virus.fa')

# q1
print len(naive('AGGT', phix_genome)) + len(naive('ACCT', phix_genome))
print len(naive_with_rc('AGGT', phix_genome))

# q2
print len(naive('TTAA', phix_genome))
print len(naive_with_rc('TTAA', phix_genome))

# q3
print naive('ACTAAGT', phix_genome)[0]
print naive('ACTTAGT', phix_genome)[0]
print naive(reverseComplement('ACTAAGT'), phix_genome)[0]

# q4
print naive('AGTCGA', phix_genome)[0]
print naive('TCGACT', phix_genome)[0]
print naive(reverseComplement('AGTCGA'), phix_genome)[0]

print naive_with_rc('AGTCGA', phix_genome)[0]

# q5
print len(naive_2mm('TTCAAGCC', phix_genome))

# q6
print naive_2mm('AGGAGGTT', phix_genome)[0]

def toQual(s):
    return ord(s) - 33

reads, qualities = readFastq('ERR037900_1.first1000.fastq')
# q7 use quality score
qlist = [0] * len(qualities[0])
for qual in qualities:
    for i in range(len(qual)):
        qlist[i] += toQual(qual[i])
for i in range(len(qlist)):
    qlist[i] = float(qlist[i]) / len(qualities)

print qlist.index(min(qlist)) #66

# q7 use base ratio
gc = [0.] * len(reads[0])
total = [0] * len(reads[0])
for read in reads:
    for i in range(len(read)):
        if read[i] == 'C' or read[i] == 'G':
            gc[i] += 1
        total[i] += 1
for i in range(len(gc)):
    gc[i] /= total[i]
    if gc[i] < 0.2:
        print i


