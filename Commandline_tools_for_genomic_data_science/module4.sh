#build bowtie2 index
bowtie2-build -f athal_chr.fa athal 1>> bowtie2.build.index.txt

#align reads with tophat
WORKDIR=/Users/Kai/Lu/MOOC/coursera/genomic_certificate/commandline_tool/module4/gencommand_proj4/
mkdir -p $WORKDIR/tophat2.Day8
tophat -o $WORKDIR/tophat2.Day8 $WORKDIR/athal $WORKDIR/Day8.fastq 2>> $WORKDIR/tophat2.Day8/tophat2.Day8.txt
mkdir -p $WORKDIR/tophat2.Day16
tophat -o $WORKDIR/tophat2.Day16 $WORKDIR/athal $WORKDIR/Day16.fastq 2>> $WORKDIR/tophat2.Day16/tophat2.Day16.txt

samtools view tophat2.Day8/accepted_hits.bam|wc -l
#q1   63845
samtools view tophat2.Day16/accepted_hits.bam|wc -l
#q2	  58398

#q3: 63489 q4: 57951 q5: 63133 q6: 57504

samtools view $WORKDIR/tophat2.Day8/accepted_hits.bam | cut -f6 | grep "N" | wc -l
#q7 8596
samtools view $WORKDIR/tophat2.Day16/accepted_hits.bam | cut -f6 | grep "N" | wc -l
#q8 10695

samtools view $WORKDIR/tophat2.Day8/unmapped.bam | wc -l
#q9 84
samtools view $WORKDIR/tophat2.Day16/unmapped.bam | wc -l
#q10 34

#assemble alignments with cufflinks
mkdir -p $WORKDIR/cufflinks.Day8
cufflinks -L Day8 -o $WORKDIR/cufflinks.Day8/ $WORKDIR/tophat2.Day8/accepted_hits.bam 2>> $WORKDIR/cufflinks.Day8/cufflinks.Day8.txt
mkdir -p $WORKDIR/cufflinks.Day16
cufflinks -L Day16 -o $WORKDIR/cufflinks.Day16/ $WORKDIR/tophat2.Day16/accepted_hits.bam 2>> $WORKDIR/cufflinks.Day16/cufflinks.Day16.txt

cut -f9 cufflinks.Day8/transcripts.gtf | cut -d ' ' -f2 | uniq | wc -l
#q11 186
cut -f9 cufflinks.Day16/transcripts.gtf | cut -d ' ' -f2 | uniq | wc -l
#q12 80

cut -f3 cufflinks.Day8/transcripts.gtf | grep 'transcript' | wc -l
#q13 192
cut -f3 cufflinks.Day16/transcripts.gtf | grep 'transcript' | wc -l
#q14 92


cut -f9 cufflinks.Day8/transcripts.gtf | cut -d ' ' -f4 | grep -E '.2";$' | uniq | wc -l
#q15 186 -6 = 180

cut -f9 cufflinks.Day16/transcripts.gtf | cut -d ' ' -f4 | grep -E '.2";$' | uniq | wc -l
#q16 80 - 11 = 69

grep "exon" cufflinks.Day8/transcripts.gtf | cut -f9 | cut -d ' ' -f6 | grep -E '^"2";$' | wc -l
#q17 192 - 73 = 119

grep "exon" cufflinks.Day16/transcripts.gtf | cut -f9 | cut -d ' ' -f6 | grep -E '^"2";$' | wc -l
#q18 92 - 68 = 24

#q19 73 q20 68

cuffcompare -r athal_genes.gtf -R -o cuffcmp.Day8 cufflinks.Day8/transcripts.gtf
cuffcompare -r athal_genes.gtf -R -o cuffcmp.Day16 cufflinks.Day16/transcripts.gtf

cut -f4 cuffcmp.Day8.tracking | grep "=" | wc -l
#16
cut -f4 cuffcmp.Day16.tracking | grep "=" | wc -l
#36

grep "AT4G20240" cuffcmp.Day8.tracking |wc -l
#2
grep "AT4G20240" cuffcmp.Day16.tracking | wc -l
#0

cut -f4 cuffcmp.Day8.tracking | grep "c" | wc -l
#133
cut -f4 cuffcmp.Day16.tracking | grep "c" | wc -l
#21

cut -f4 cuffcmp.Day8.tracking | grep "j" | wc -l
#14
cut -f4 cuffcmp.Day16.tracking | grep "j" | wc -l
#22

cut -f4 cuffcmp.Day8.tracking | grep "i" | wc -l
#4
cut -f4 cuffcmp.Day16.tracking | grep "i" | wc -l
#1

cuffmerge -g athal_genes.gtf assembly_GTF_list.txt 

cut -f9 merged_asm/merged.gtf | cut -d ' ' -f2 | uniq | wc -l
#129
cut -f9 merged_asm/merged.gtf | cut -d ' ' -f4 | uniq | wc -l
#200

cuffdiff --max-bundle-frags 500000 -o $WORKDIR/cuffdiff merged_asm/merged.gtf tophat2.Day8/accepted_hits.bam tophat2.Day16/accepted_hits.bam

cut -f2 cuffdiff/gene_exp.diff | uniq | wc -l
#130-1=129 (subtract the header)

less cuffdiff/gene_exp.diff | grep yes  | wc -l
#4

less cuffdiff/isoform_exp.diff | grep yes  | wc -l
#5






