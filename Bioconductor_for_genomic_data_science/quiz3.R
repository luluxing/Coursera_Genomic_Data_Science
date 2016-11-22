library(ALL)
data(ALL)
expression.data = exprs(ALL)
mean(expression.data[, 5]) 
#q1 5.629627

library(biomaRt)
mart = useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
hg19.ensembl = useDataset("hsapiens_gene_ensembl", mart)
vals = rownames(expression.data)
atts = c("ensembl_gene_id", "affy_hg_u95av2")
ensmbl_id.feat = getBM(attributes = atts, filters = "affy_hg_u95av2", values = vals, mart = hg19.ensembl)
length(which(table(ensmbl_id.feat[,2])>1))
#q2 1045

total.attri = listAttributes(hg19.ensembl)
head(total.attri)
ensmbl_id.feat2 = getBM(attributes = c('chromosome_name', "ensembl_gene_id", "affy_hg_u95av2"), filters = "affy_hg_u95av2", values = vals, mart = hg19.ensembl)
autosomes = as.character(1:22)
ensmbl_id.feat3 = ensmbl_id.feat2[ensmbl_id.feat2[,1] %in% autosomes, ]
length(table(ensmbl_id.feat3[,3]))
#q3 11016

library(minfiData)
data(MsetEx)
methyl = getMeth(MsetEx)
mean(methyl[, '5723646052_R04C01'])
#q4 7228.277

library(GEOquery)
gse = getGEO(filename = 'GSE788_family.soft.gz')
mean(Table(GSMList(gse)[[2]])$VALUE)
#q5 756.432

library(airway)
data(airway)
mean(colData(airway)$avgLength)
#q6 113.75

exprs_data = assay(airway, "counts")
sum(exprs_data[,'SRR1039512']>=1)
#q7 25699

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exon.txdb = exons(txdb)
autosomes = paste0("chr", c(1:22))
exon.txdb = keepSeqlevels(exon.txdb, autosomes)
exon.txdb.rename = renameSeqlevels(exon.txdb, mapSeqlevels(autosomes, 'NCBI'))
airway_txdb.overlap = subsetByOverlaps(airway, exon.txdb.rename)
#q8 2627

total_reads = sum(exprs_data[,'SRR1039508'])
airway_txdb.overlap_data = assay(airway_txdb.overlap, "counts")
overlap_reads = sum(airway_txdb.overlap_data[,'SRR1039508'])
overlap_reads/total_reads
#q9 0.9004193

library(AnnotationHub)
ah = AnnotationHub()
query(ah, c('E096','H3K4me3'))
ah_que = ah[['AH30596']]
ah_que = keepSeqlevels(ah_que, autosomes)
ah_que.rename = renameSeqlevels(ah_que, mapSeqlevels(autosomes, 'NCBI'))

auto_ncbi = extractSeqlevelsByGroup(species="Homo sapiens", style="NCBI", group="auto")
airway_txdb.overlap = keepSeqlevels(airway_txdb.overlap, auto_ncbi)

airway.overlap.prom = promoters(airway_txdb.overlap)
airway_ah.overlap = subsetByOverlaps(airway.overlap.prom, ah_que.rename)
airway_ah.overlap_data = assay(airway_ah.overlap, "counts")
median(airway_ah.overlap_data[,'SRR1039508'])



