library(Biobase)



#q1
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata2 = log2(edata + 1)
edata3 = edata2 - rowMeans(edata2)
edata1 = edata
svd1 = svd(edata1)
svd2 = svd(edata2)
svd3 = svd(edata3)
(svd1$d^2/sum(svd1$d^2))[1] # 0.8873421
(svd2$d^2/sum(svd2$d^2))[1] # 0.9737781
(svd3$d^2/sum(svd3$d^2))[1] # 0.3463729

#q2
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata2 = log2(edata + 1)
edata3 = edata2 - rowMeans(edata2)
edata1 = edata
svd1 = svd(edata1)
svd2 = svd(edata2)
svd3 = svd(edata3)
set.seed(333)
edata3.kms = kmeans(t(edata3), 2)
cor(svd3$v[,1], as.numeric(edata3.kms$cluster)) # -0.8678247

#q3
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
lm1 = lm(edata[1,] ~ pdata_bm$num.tech.reps)
plot(pdata_bm$num.tech.reps, edata[1,])
abline(lm1)

#q4
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
lm2 = lm(edata[1,] ~ pdata_bm$gender + pdata_bm$age)
library(broom)
tidy(lm2) # -23.9133

#q5
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(as.matrix(edata) + 1)
mod1 = model.matrix( ~ pdata$population)
fit1 = lm.fit(mod1, t(edata))
dim(fit1$residuals)
# 129 52580
dim(fit1$effects)
# 129 52580
dim(fit1$coefficients)
# 2 52580

#q7
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
na.cols = pdata_bm$sample.id[is.na(pdata_bm$age == 'NA')]
edata.rm.age.na = edata[, !(colnames(edata) %in% na.cols)]
pdata_bm.rm.na = pdata_bm[!(rownames(pdata_bm) %in% na.cols),]
mod2 = model.matrix( ~ pdata_bm.rm.na$age)
fit.limma = lmFit(edata.rm.age.na, mod2)
fit.limma$coefficients[1000,]
# (Intercept) pdata_bm$age 
#  2469.87375    -27.61178 
plot(pdata_bm.rm.na$age, edata.rm.age.na[1000,])

#q8
mod2_adj = model.matrix( ~ pdata_bm.rm.na$age + as.factor(pdata_bm.rm.na$tissue.type))
fit.limma.adj = lmFit(edata.rm.age.na, mod2_adj)

#q10
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
set.seed(33353)
library(sva)
library(snpStats)
edata.trans = log2(edata + 1)
edata.trans = edata.trans[rowMeans(edata.trans) >= 1, ]
na.cols = pdata_bm$sample.id[is.na(pdata_bm$age == 'NA')]
edata.rm.na = edata.trans[, !(colnames(edata.trans) %in% na.cols)]
pdata_bm.rm.na = pdata_bm[!(rownames(pdata_bm) %in% na.cols),]

mod = model.matrix( ~ pdata_bm.rm.na$age, data = pdata_bm.rm.na)
mod0 = model.matrix(~1, data=pdata_bm.rm.na)
sva1 = sva(edata.rm.na, mod, mod0, n.sv=1)
cor(sva1$sv, pdata_bm.rm.na$age)
