library(AnnotationHub)
ah = AnnotationHub()
ah_hg19 = subset(ah, species == 'Homo sapiens')
ah_hg19_cpg = query(ah_hg19, c('CpG Islands'))
gr_hg19_cpg = ah_hg19_cpg[[1]]
autosomes = paste ("chr", 1:22, sep = '')
seqlevels(gr_hg19_cpg, force = TRUE) = autosomes
gr_hg19_cpg # 26641

gr_hg19_cpg_chr4 = gr_hg19_cpg
seqlevels(gr_hg19_cpg_chr4, force = TRUE) = 'chr4'
gr_hg19_cpg_chr4 # 1031

qhs_h3k4m3 = query(ah, c("H3K4me3", "E003"))
qhs_h3k4m3 = qhs_h3k4m3[['AH29884']]
seqlevels(qhs_h3k4m3, force = TRUE) = autosomes
sum(width(reduce(qhs_h3k4m3))) # 41135164

qhs_h3k27m3 = query(ah, c("H3K27me3", "E003"))
qhs_h3k27m3 = qhs_h3k27m3[['AH29892']]
seqlevels(qhs_h3k27m3, force = TRUE) = autosomes
mean(qhs_h3k27m3$signalValue) # 4.770728

sum(width(intersect(qhs_h3k4m3, qhs_h3k27m3))) # 10289096

bivalent = intersect(qhs_h3k4m3, qhs_h3k27m3)
length(subsetByOverlaps(bivalent, gr_hg19_cpg, ignore.strand=TRUE)) / length(bivalent) # 0.5383644

sum(width(intersect(bivalent, gr_hg19_cpg))) / sum(width(reduce(gr_hg19_cpg))) # 0.241688

CpG_10k <- resize(unlist(gr_hg19_cpg), width = 20000 + width(unlist(gr_hg19_cpg)), fix = "center")
sum(width(intersect(bivalent, CpG_10k))) # 9782086

qhs_hg19 = ah_hg19[['AH5018']]
seqlevels(qhs_hg19, force = TRUE) = autosomes
genome_size = sum(as.numeric(seqlengths(qhs_hg19)))
sum(width(reduce(gr_hg19_cpg))) / genome_size # 0.007047481

inout = matrix(0, ncol = 2, nrow = 2)
colnames(inout) = c('in', 'out')
rownames(inout) = c('in', 'out')
inout[1,1] = sum(width(intersect(bivalent, gr_hg19_cpg, ignore.strand=TRUE)))
inout[1,2] = sum(width(setdiff(bivalent, gr_hg19_cpg, ignore.strand=TRUE)))
inout[2,1] = sum(width(setdiff(gr_hg19_cpg, bivalent, ignore.strand=TRUE)))
inout[2,2] = genome_size - sum(inout)
oddsratio = inout[1,1] * inout[2,2] / (inout[2,1] * inout[1,2])
oddsratio # 169.0962