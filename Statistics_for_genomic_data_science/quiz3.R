library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA
glm1 = glm(status ~ snp3, family= "binomial")
#q1
# > tidy(glm1)
#          term   estimate std.error  statistic   p.value
# 1 (Intercept)  0.1772078 0.2199584  0.8056424 0.4204491
# 2        snp3 -0.1579378 0.1877400 -0.8412582 0.4002033
# > tidy(lm1)
#          term    estimate  std.error  statistic      p.value
# 1 (Intercept)  0.54421636 0.05489222  9.9142712 3.749973e-22
# 2        snp3 -0.03940072 0.04681654 -0.8415982 4.002165e-01

snp10 = as.numeric(snpdata[,10])
snp10[snp10==0] = NA
snp10_ress = (snp10 == 2)
glm10_ress = glm(status ~ snp10_ress, family="binomial")
glm10_add = glm(status ~ snp10, family= "binomial")
glm10_ress$fitted.values
glm10_add$fitted.values
#q3
#No, in all cases, the fitted values are near 0.5 and there 
#are about an equal number of cases and controls in each group. 
#This is true regardless of whether you fit a recessive or additive model.

effects = rep(0,2851)
for (i in 1:ncol(snpdata)) {
	tmp = as.numeric(snpdata[,i])
	tmp[tmp==0] = NA
	glm = glm(status ~ tmp, family="binomial")
	effects[i] = tidy(glm)$statistic[2]
}
#q4
# > mean(effects)
# [1] 0.007155377
# > min(effects)
# [1] -4.251469
# > max(effects)
# [1] 3.900891

coeffi = rep(0,2851)
for (i in 1:ncol(snpdata)) {
	tmp = as.numeric(snpdata[,i])
	tmp[tmp==0] = NA
	glm = glm(status ~ tmp, family="binomial")
	coeffi[i] = glm$coefficients
}
squared_coeff = coeffi^2
glm_all = snp.rhs.tests(status ~ 1, snp.data=sub.10)
chi.glm_all = chi.squared(glm_all)
#q5
# cor(squared_coeff, chi.glm_all)
# [1] 0.04566279
# > 0.99. They are both testing for the same association 
# using the same additive regression model on the logistic 
# scale but using slightly different tests.

library(genefilter)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(as.matrix(edata) + 1)
tstatic = rowttests(edata, pdata$population)
fstatic = rowFtests(edata, pdata$population)
#q6
# head(tstatic$p.value)
# head(fstatic$p.value)
# head(tstatic$statistic)
# head(fstatic$statistic)
# You get the same p-value but different statistics. 
# This is because the F-statistic and t-statistic 
# test the same thing when doing a two group test and 
# one is a transform of the other.

library(DESeq2)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
de = DESeqDataSetFromMatrix(edata, pdata, ~population)
glm_all_nb = DESeq(de)
deseq_result = results(glm_all_nb)

edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~pdata$population)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)
ebayes_limma.topt = topTable(ebayes_limma,sort="none",n=Inf)
cor(deseq_result$stat, ebayes_limma.topt$t)
#q7 0.9278701
MA = new("MAList")
MA$M = deseq_result$stat
MA$A = ebayes_limma.topt$t
limma::plotMA(MA)
#There are more differences for the small statistics.

despadj = p.adjust(deseq_result$pvalue, method = 'BH')
limpadj = p.adjust(ebayes_limma.topt$P.Value, method = 'BH')
#q8
# > sum(despadj<=0.05)
# [1] 1995
# > sum(limpadj<=0.05)
# [1] 2807

