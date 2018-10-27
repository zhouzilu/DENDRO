#!/bin/bash

# 1. Index the genome template. Use STAR 2 pass. STAR manuel http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf
# This step only needs to be run once
# pathtoSTAR/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 16 --genomeDir pathtoSTAR/STAR-2.5.3a/star_genome --genomeFastaFiles pathtoSTAR/STAR_hg19/ucsc.hg19.fasta

SAMPLE=$1
echo "${SAMPLE}"

# 2. Mapping by STAR 2 pass
# 2.1 first pass
mv ../fastq/${SAMPLE}_1.fastq ../fastq/${SAMPLE}_2.fastq .
pathtoSTAR/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --genomeDir pathtoSTAR/STAR-2.5.3a/star_genome --runThreadN 24 --readFilesIn ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq --outFileNamePrefix ${SAMPLE}
# 2.2 combine all the SJ.out.tab to regenerate reference. You need to list all the SJ.out.tab files, which is painful.
pathtoSTAR/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --genomeDir pathtoSTAR/STAR-2.5.3a/star_genome --genomeFastaFiles pathtoSTAR/STAR_hg19/ucsc.hg19.fasta --sjdbFileChrStartEnd SRR2973275SJ.out.tab SRR2973365SJ.out.tab SRR2973380SJ.out.tab SRR5023592SJ.out.tab SRR5023607SJ.out.tab SRR5023622SJ.out.tab SRR5023637SJ.out.tab SRR2973351SJ.out.tab SRR2973366SJ.out.tab SRR2973381SJ.out.tab SRR5023593SJ.out.tab SRR5023608SJ.out.tab SRR5023623SJ.out.tab SRR5023638SJ.out.tab SRR2973352SJ.out.tab SRR2973367SJ.out.tab SRR2973382SJ.out.tab SRR5023594SJ.out.tab SRR5023609SJ.out.tab SRR5023624SJ.out.tab SRR5023639SJ.out.tab SRR2973353SJ.out.tab SRR2973368SJ.out.tab SRR2973383SJ.out.tab SRR5023595SJ.out.tab SRR5023610SJ.out.tab SRR5023625SJ.out.tab SRR5023640SJ.out.tab SRR2973354SJ.out.tab SRR2973369SJ.out.tab SRR5023378SJ.out.tab SRR5023596SJ.out.tab SRR5023611SJ.out.tab SRR5023626SJ.out.tab SRR5023641SJ.out.tab SRR2973355SJ.out.tab SRR2973370SJ.out.tab SRR5023442SJ.out.tab SRR5023597SJ.out.tab SRR5023612SJ.out.tab SRR5023627SJ.out.tab SRR5023642SJ.out.tab SRR2973356SJ.out.tab SRR2973371SJ.out.tab SRR5023443SJ.out.tab SRR5023598SJ.out.tab SRR5023613SJ.out.tab SRR5023628SJ.out.tab SRR5023643SJ.out.tab SRR2973357SJ.out.tab SRR2973372SJ.out.tab SRR5023444SJ.out.tab SRR5023599SJ.out.tab SRR5023614SJ.out.tab SRR5023629SJ.out.tab SRR5023644SJ.out.tab SRR2973358SJ.out.tab SRR2973373SJ.out.tab SRR5023445SJ.out.tab SRR5023600SJ.out.tab SRR5023615SJ.out.tab SRR5023630SJ.out.tab SRR5023645SJ.out.tab SRR2973359SJ.out.tab SRR2973374SJ.out.tab SRR5023586SJ.out.tab SRR5023601SJ.out.tab SRR5023616SJ.out.tab SRR5023631SJ.out.tab SRR2973360SJ.out.tab SRR2973375SJ.out.tab SRR5023587SJ.out.tab SRR5023602SJ.out.tab SRR5023617SJ.out.tab SRR5023632SJ.out.tab SRR2973361SJ.out.tab SRR2973376SJ.out.tab SRR5023588SJ.out.tab SRR5023603SJ.out.tab SRR5023618SJ.out.tab SRR5023633SJ.out.tab SRR2973362SJ.out.tab SRR2973377SJ.out.tab SRR5023589SJ.out.tab SRR5023604SJ.out.tab SRR5023619SJ.out.tab SRR5023634SJ.out.tab SRR2973363SJ.out.tab SRR2973378SJ.out.tab SRR5023590SJ.out.tab SRR5023605SJ.out.tab SRR5023620SJ.out.tab SRR5023635SJ.out.tab SRR2973364SJ.out.tab SRR2973379SJ.out.tab SRR5023591SJ.out.tab SRR5023606SJ.out.tab SRR5023621SJ.out.tab SRR5023636SJ.out.tab --sjdbOverhang 75 --runThreadN 24
# 2.3 second pass
pathtoSTAR/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --genomeDir pathtoSTAR/STAR-2.5.3a/star_genome --runThreadN 24 --readFilesIn ../fastq/${SAMPLE}_1.fastq ../fastq/${SAMPLE}_2.fastq --outFileNamePrefix ${SAMPLE}

# 3. Convert .sam to .bam and sort
echo "$1 sam to bam and sort..."
pathtosamtools/samtools-1.5/samtools view -bS ${SAMPLE}Aligned.out.sam > ${SAMPLE}.bam
java -jar pathtopicard/picard/picard.jar SortSam INPUT=${SAMPLE}.bam OUTPUT=${SAMPLE}.sorted.bam SORT_ORDER=coordinate
# rm ${SAMPLE}.bam

# 4. Add read group
echo "$1 add read group..."
java -jar pathtopicard/picard/picard.jar AddOrReplaceReadGroups INPUT=${SAMPLE}.sorted.bam OUTPUT=${SAMPLE}.sorted.rg.bam RGID=${SAMPLE} RGLB=trancriptome RGPL=ILLUMINA RGPU=machine RGSM=${SAMPLE}
pathtosamtools/samtools-1.5/samtools index ${SAMPLE}.sorted.rg.bam
rm ${SAMPLE}.sorted.bam

# 5. Dedup
echo "$1 dedup..."
java -jar pathtopicard/picard/picard.jar MarkDuplicates INPUT=${SAMPLE}.sorted.rg.bam OUTPUT=${SAMPLE}.sorted.rg.dedup.bam METRICS_FILE=${SAMPLE}.sorted.rg.dedup.metrics.txt PROGRAM_RECORD_ID= MarkDuplicates PROGRAM_GROUP_VERSION=null PROGRAM_GROUP_NAME=MarkDuplicates
java -jar pathtopicard/picard/picard.jar BuildBamIndex INPUT=${SAMPLE}.sorted.rg.dedup.bam
rm ${SAMPLE}.sorted.rg.bam

# 6. SplitNCigarReads
java -jar pathtogatk/gatk/GenomeAnalysisTK.jar -T SplitNCigarReads -R pathtoSTAR/STAR_hg19/ucsc.hg19.fasta -I ${SAMPLE}.sorted.rg.dedup.bam -o ${SAMPLE}.sorted.rg.dedup.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
rm ${SAMPLE}.sorted.rg.dedup.bam

# 7. Realign
echo "$1 realign..."
java -jar pathtogatk/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R pathtoSTAR/STAR_hg19/ucsc.hg19.fasta -I ${SAMPLE}.sorted.rg.dedup.split.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o ${SAMPLE}.sorted.rg.dedup.split.target_intervals.list
java -jar pathtogatk/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R pathtoSTAR/STAR_hg19/ucsc.hg19.fasta -I ${SAMPLE}.sorted.rg.dedup.split.bam -targetIntervals ${SAMPLE}.sorted.rg.dedup.split.target_intervals.list -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o ${SAMPLE}.sorted.rg.dedup.split.realigned.bam
rm ${SAMPLE}.sorted.rg.dedup.split.bam


# 8. Recalibrate
echo "$1 recalibrate..."
java -jar pathtogatk/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R pathtoSTAR/STAR_hg19/ucsc.hg19.fasta -I ${SAMPLE}.sorted.rg.dedup.split.realigned.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o ${SAMPLE}.sorted.rg.dedup.split.realigned.recal_data.table
java -jar pathtogatk/gatk/GenomeAnalysisTK.jar -T PrintReads -R pathtoSTAR/STAR_hg19/ucsc.hg19.fasta -I ${SAMPLE}.sorted.rg.dedup.split.realigned.bam -BQSR ${SAMPLE}.sorted.rg.dedup.split.realigned.recal_data.table -o ${SAMPLE}.sorted.rg.dedup.split.realigned.recal.bam
pathtosamtools/samtools-1.5/samtools index ${SAMPLE}.sorted.rg.dedup.split.realigned.recal.bam
rm ${SAMPLE}.sorted.rg.dedup.split.target_intervals.list ${SAMPLE}.sorted.rg.dedup.split.realigned.bam ${SAMPLE}.sorted.rg.dedup.metrics.txt

# 9. GATK HaplotypeCaller
echo "$1 variants calling... map"
java -jar pathtogatk/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller -R pathtoSTAR/STAR_hg19/ucsc.hg19.fasta -I ${SAMPLE}.sorted.rg.dedup.split.realigned.recal.bam -dontUseSoftClippedBases -ERC GVCF -stand_call_conf 20 -o GVCF/${SAMPLE}.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf

# 10. GATK HaplotypeCaller reduce GVCFs
# You also need specify each g.vcf files......
java -jar pathtogatk/gatk/GenomeAnalysisTK.jar -T GenotypeGVCFs -R pathtostar/STAR_hg19/ucsc.hg19.fasta \
-V SRR2973275.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973351.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973352.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973353.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973354.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973355.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973356.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973357.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973358.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973359.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973360.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973361.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973362.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973363.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973364.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973365.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973366.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973367.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973368.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973369.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973371.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973372.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973373.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973374.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973375.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973376.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973377.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973378.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973379.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973380.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973381.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973382.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973383.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023378.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023442.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023587.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023588.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023589.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023590.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023591.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023592.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023593.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023595.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023596.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023597.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023598.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023599.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023600.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023601.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023602.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023603.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023604.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023605.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023607.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023608.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023609.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023610.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023611.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023612.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023613.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023614.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023615.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023616.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023617.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023619.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023620.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023621.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023622.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023623.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023624.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023626.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023627.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023628.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023629.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023630.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023632.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023633.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023635.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023636.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023637.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023640.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023641.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023642.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023643.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023644.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR5023645.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973276.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973384.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973385.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973386.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973387.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973388.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973389.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973390.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973391.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973392.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973393.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973394.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973395.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973396.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973397.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973398.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973399.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973400.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973401.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973402.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973403.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973404.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973405.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973406.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973407.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973408.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973409.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973410.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973411.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973412.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973413.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973414.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973416.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973417.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973418.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973420.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973421.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973422.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973423.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973424.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973426.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973427.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973429.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973430.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973431.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973433.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973434.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973435.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-V SRR2973436.sorted.rg.dedup.realigned.recal.raw.snps.indels.cof20.erc.g.vcf \
-o BC030903LN.sorted.rg.dedup.split.realigned.recal.bam.cof20.GVCFs.jointcalls.g.vcf
