grep '>' wu_0.v7.fas | wc -l # 7

bowtie2-build wu_0.v7.fas wu_0

wc -l wu_0_A_wgs.fastq # 147354

# bowtie2 -x wu_0 wu_0_A_wgs.fastq -S wu_0.bt2.sam ## different version leads to different results
bowtie2 -x wu_0 wu_0_A_wgs.fastq -S wu_0.bt2.oldversion.sam
samtools view -b wu_0.bt2.oldversion.sam | samtools sort - > wu_0.bt2.sorted.oldversion.bam
samtools mpileup -v -u -f wu_0.v7.fas wu_0.bt2.sorted.oldversion.bam > wu_0.bt2.oldversion.vcf
bcftools view -H wu_0.bt2.oldversion.vcf | cut -f4 | grep '^A$' | wc -l # 1150985


bowtie2 --local -x wu_0 wu_0_A_wgs.fastq -S wu_0.bt2.local.sam

# samtools view wu_0.bt2.sam | cut -f3 | grep -v '*' | wc -l # 137718
samtools view wu_0.bt2.oldversion.sam | cut -f3 | grep -v '*' | wc -l # 137719


samtools view wu_0.bt2.local.sam | cut -f3 | grep -v '*' | wc -l #  141039

samtools view wu_0.bt2.sam | cut -f6 | grep -E 'I|D'| wc -l # 2782

# samtools view wu_0.bt2.local.sam | cut -f6 | grep -E 'I|D'| wc -l # 2612
bowtie2 --local -x wu_0 wu_0_A_wgs.fastq -S wu_0.bt2.local.oldversion.sam
samtools view wu_0.bt2.local.oldversion.sam | cut -f6 | grep -E 'I|D'| wc -l
#2613

samtools view -b wu_0.bt2.sam | samtools sort - > wu_0.bt2.sorted.bam

samtools index wu_0.bt2.sorted.bam

samtools mpileup -v -u -f wu_0.v7.fas wu_0.bt2.sorted.bam > wu_0.bt2.vcf

bcftools view -H wu_0.bt2.vcf | cut -f1 | grep 'Chr3' | wc -l # 360295

# bcftools view -H wu_0.bt2.vcf | cut -f4 | grep '^A$' | wc -l # 1150980

bcftools view -H wu_0.bt2.vcf | cut -f8 | grep 'DP=20' | wc -l # 1816

bcftools view -H wu_0.bt2.vcf | cut -f8 | grep 'INDEL' | wc -l # 1972

bcftools view -H wu_0.bt2.vcf | cut -f1,2 | grep 'Chr1'| grep '175672' | wc -l # 2

samtools mpileup -g -f wu_0.v7.fas wu_0.bt2.sorted.bam > wu_0.bt2.bcf

bcftools call -v -m -O v -o wu_0.bt2.variantsonlly.vcf wu_0.bt2.bcf

bcftools view -H wu_0.bt2.variantsonlly.vcf | cut -f1 | grep 'Chr3' | wc -l # 398

bcftools view -H wu_0.bt2.variantsonlly.vcf | cut -f4,5 | grep -E '^A\tT' | wc -l # 392

bcftools view -H wu_0.bt2.variantsonlly.vcf | cut -f8 | grep 'INDEL' | wc -l # 320

bcftools view -H wu_0.bt2.variantsonlly.vcf | cut -f8 | grep 'DP=20' | wc -l # 2

bcftools view -H wu_0.bt2.variantsonlly.vcf | grep 'Chr3'| grep '11937923' | less # SNP
