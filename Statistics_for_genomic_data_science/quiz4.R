library(goseq)
supportedGenomes()
#q1 UCSC mm9

library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
edata = log2(edata+1)

mod = model.matrix(~pdata_bot$strain)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)
ebayes_limma.topt = topTable(ebayes_limma,sort="none",n=Inf)
limpadj = p.adjust(ebayes_limma.topt$P.Value, method = 'BH')
#q2
# > sum(limpadj<=0.05)
# [1] 223
# featureNames(fdata_bot[ebayes_limma.topt$adj.P.Val<=0.05])[1]
# [1] "ENSMUSG00000000402"

# genes = as.integer(p.adjust(ebayes_limma.topt$P.Value[ebayes_limma.topt$logFC!=0], method="BH")<.05)
genes = as.integer(p.adjust(ebayes_limma.topt$P.Value, method="BH")<.05)
names(genes)=row.names(ebayes_limma.topt)
pwf = nullp(genes,"mm9","ensGene")
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
pvals = goseq(pwf, 'mm9','ensGene')
head(pvals)
#q3 GO:0004888
#q4 transmembrane signaling receptor activity

library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]

mod = model.matrix(~pdata_bot$strain)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)
ebayes_limma.topt = topTable(ebayes_limma,sort="none",n=Inf)
genes = as.integer(p.adjust(ebayes_limma.topt$P.Value, method="BH")<=.05)
names(genes)=row.names(ebayes_limma.topt)
pwf = nullp(genes,"mm9","ensGene")
pvals = goseq(pwf, 'mm9','ensGene')

mod2 = model.matrix(~ pdata_bot$strain + as.factor(pdata_bot$lane.number))
fit_limma2 = lmFit(edata, mod2)
ebayes_limma2 = eBayes(fit_limma2)
ebayes_limma2.topt = topTable(ebayes_limma2,sort="none",n=Inf)
genes2 = as.integer(p.adjust(ebayes_limma2.topt$P.Value, method="BH")<=.05)
names(genes2)=row.names(ebayes_limma2.topt)
pwf2 = nullp(genes2,"mm9","ensGene")
pvals2 = goseq(pwf2, 'mm9','ensGene')

unadjgo = pvals$category[1:10]
adjgo = pvals2$category[1:10]
intersect(adjgo, unadjgo)
#q5 2




