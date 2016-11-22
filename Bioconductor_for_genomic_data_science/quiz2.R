library(BSgenome)
library(Biostrings)
library(GenomicFeatures)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("BSgenome.Hsapiens.UCSC.hg19")

#q1
letterFrequency(Hsapiens$chr22, 'GC') / letterFrequency(Hsapiens$chr22, 'ATCG') # 0.4798807

#q2
library(AnnotationHub)
ah = AnnotationHub()
qhs_h3k27m3 = query(ah, c("H3K27me3", "E003"))
qhs_h3k27m3 = qhs_h3k27m3[['AH29892']]
seqlevels(qhs_h3k27m3, force = TRUE) = 'chr22'
hs_chr22 = Views(Hsapiens, qhs_h3k27m3)
sum(unlist(letterFrequency(hs_chr22, 'GC', as.prob = TRUE))) / length(letterFrequency(hs_chr22, 'GC', as.prob = TRUE))
# 0.528866

#q3
cor(letterFrequency(hs_chr22, 'GC', as.prob = TRUE), qhs_h3k27m3$signalValue) # 0.004467924

#q4
fc_h3k27m3 = query(ah, c("H3K27me3", "E003", "fc.signal"))
fc_h3k27m3 = fc_h3k27m3[["AH32033"]]
library(rtracklayer)
gr.rle = import(fc_h3k27m3, which=GRanges('chr22', ranges=IRanges(1, 10^8)), as = 'Rle')
gr.chr22.rle = gr.rle$chr22
fc.chr22 = Views(gr.chr22.rle, start=start(qhs_h3k27m3), end=end(qhs_h3k27m3))
fc.chr22.mean = mean(fc.chr22)
cor(fc.chr22.mean, qhs_h3k27m3$signalValue) # 0.9149614

#q5
sum(gr.chr22.rle>1) # 10914671

#q6
stemcell_h3k27m3 = query(ah, c("H3K27me3", "E055", "fc.signal"))
stemcell_h3k27m3 = stemcell_h3k27m3[["AH32470"]]
stemcell_h3k27m3.rle = import(stemcell_h3k27m3, which=GRanges('chr22', ranges=IRanges(1, 10^8)), as = 'Rle')
stemcell_h3k27m3.chr22.rle = stemcell_h3k27m3.rle$chr22
stemcell.ir = as()
E003.lower0.5 = slice(gr.chr22.rle, upper=0.5, includeUpper=TRUE)
E055.higher2 = slice(stemcell_h3k27m3.chr22.rle, lower=2, includeLower=TRUE)
E003.lower0.5.ir = as(E003.lower0.5, 'IRanges')
E055.higher2.ir = as(E055.higher2, 'IRanges')
sum(width(intersect(E003.lower0.5.ir, E055.higher2.ir))) # 1869937

#q7
hg19.cpg = query(subset(ah, species == 'Homo sapiens'), c('CpG Islands'))
hg19.cpg = hg19.cpg[["AH5086"]]
hg19.chr22.cpg = subset(hg19.cpg, seqnames == "chr22")
hg19.chr22.cpg.vi = Views(Hsapiens, hg19.chr22.cpg)
freq.C = letterFrequency(hg19.chr22.cpg.vi, 'C')
freq.G = letterFrequency(hg19.chr22.cpg.vi, 'G')
obse.CG = dinucleotideFrequency(hg19.chr22.cpg.vi)[,'CG']
freq.exp.CG = freq.C[,'C'] * freq.G[,'G'] / width(hg19.chr22.cpg.vi)
mean(obse.CG / freq.exp.CG) # 0.8340929

#q8
countPattern('TATAAA', Hsapiens$chr22) + countPattern('TATAAA', reverseComplement(Hsapiens$chr22))
# 27263

#q9
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
gr = GRanges(seqnames='chr22', ranges=IRanges(1, 10^8))
gr.transci.chr22 = subsetByOverlaps(transcripts(txdb), gr, ignore.strand=TRUE)
# gr.transci.chr22 = subset(transcripts(txdb), seqnames=='chr22')
proms = promoters(gr.transci.chr22, upstream=900, downstream=100)
cds = subsetByOverlaps(genes(txdb), gr, ignore.strand=TRUE)
# cds = subset(genes(txdb), seqnames=='chr22')
proms.cds = findOverlaps(proms, cds)
# cds.proms = findOverlaps(cds, proms)
# h.cds.proms = unique(queryHits(cds.proms))
h = unique(queryHits(proms.cds))
count = 0
for (i in h) {
	count = count + vcountPattern('TATAAA', DNAStringSet(Views(Hsapiens, proms[i])))
}
count

#q10
tl.chr22 = transcriptLengths(txdb, with.cds_len = TRUE)
tl.chr22 = tl.chr22[tl.chr22$cds_len > 0,]
trans.eval = proms[mcols(proms)$tx_id %in% tl.chr22$tx_id]
sum(coverage(trans.eval) > 1) # 306920





