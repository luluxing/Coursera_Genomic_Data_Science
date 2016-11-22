library(yeastRNASeq)
library(ShortRead)
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
reads = readFastq(fastqFilePath)
sum(subseq(sread(reads), start = 5, end = 5) == DNAString('A')) / length(reads)
#q1 0.363841

sum(as(quality(reads), "matrix")[, 5]) / length(reads)
#q2 28.93346

library(Rsamtools)
library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")
yeast_bamfile = BamFile(bamFilePath)
gr = GRanges(seqnames="Scchr13", ranges=IRanges(start = 800000, end = 801000))
params = ScanBamParam(which = gr, what = scanBamWhat())
yeast_aln = scanBam(yeast_bamfile, param = params)
sum(table(yeast_aln[[1]]$pos))-sum(table(yeast_aln[[1]]$pos)==1)
#q3 129

bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
yeast_bam_views = BamViews(bpaths)
gr2 = GRanges(seqnames="Scchr13", ranges=IRanges(start = 807762, end = 808068))
bamRanges(yeast_bam_views) = gr2
yeast_multi_aln = scanBam(yeast_bam_views)
l=0
for (i in 1:8) {
	l = l + length(yeast_multi_aln[[i]][[1]]$seq)
}
l/8
#q4 90.25

library(oligo)
library(GEOquery)
getGEOSuppFiles("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL")
celfiles = list.files("GSE38792/CEL", full = TRUE)
rawData = read.celfiles(celfiles)
getClass("GeneFeatureSet")
filename = sampleNames(rawData)
pData(rawData)$filename = filename
sampleNames = sub(".*_", "", filename)
sampleNames = sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) = sampleNames
pData(rawData)$group = ifelse(grepl("^OSA", sampleNames(rawData)),"OSA","Control")
normData = rma(rawData)
mean(exprs(normData)["8149273",1:8])
#q5 7.02183

library(limma)
ourData = exprs(normData)
design = model.matrix(~ pData(rawData)$group)
fit = lmFit(ourData, design)
fit = eBayes(fit)
topt = topTable(fit)
topt[which.min(topt$P.Value), 'logFC']
#q6 -0.7126484

deg_table = topTable(fit, number = length(ourData))
sum(deg_table$adj.P.Val<0.05)
#q7 0

library(minfi)
library(minfiData)
rgset_prepro = preprocessFunnorm(RGsetEx)
beta_val = assay(rgset_prepro, "Beta")
open_sea_loci = getIslandStatus(rgset_prepro)
beta_val_opensea = beta_val[open_sea_loci == "OpenSea",]
beta_means = colMeans(beta_val_opensea)
mean(beta_means[1],beta_means[2],beta_means[5])-mean(beta_means[3],beta_means[4],beta_means[6])
#q8 0.08767477

library(AnnotationHub)
ah = AnnotationHub()
ah_caco2 = query(ah, c('Caco2', 'AWG'))
DNase_hyper_sens = ah_caco2[['AH22442']]
gr_450k = granges(rgset_prepro)
subsetByOverlaps(DNase_hyper_sens, gr_450k)
#q9 40151

library(zebrafishRNASeq)
library(DESeq2)
data("zfGenes")
zfrowname = rownames(zfGenes)
zf_no_spikein = zfGenes[!startsWith(zfrowname, 'ERCC'),]
colData <- DataFrame(Treatment=c(rep("Control",3), rep("Treated", 3)), row.names=c("Ctl1","Ctl3","Ctl5","Trt9","Trt11","Trt13"))
counts = as.matrix(zf_no_spikein)
zf_assay = SummarizedExperiment(assays=list(counts=counts),colData=colData)

dds1 = DESeqDataSetFromMatrix(countData=counts, colData=colData, design = ~ Treatment)

dds = DESeqDataSet(zf_assay, design = ~ Treatment)
dds = DESeq(dds)
res = results(dds)
sum(res$padj <= 0.05)
#q10 72