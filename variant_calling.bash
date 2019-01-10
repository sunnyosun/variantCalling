#!/bin/bash                                          
#$ -S /bin/bash
#$ -cwd

module load bwa
module load picard-tools/2.6.0
module load samtools
module load gatk/3.6
module load bedtools

ref=/ifs/data/proteomics/projects/Sunny/jingchuan/yeast_genome/new/BY4741_genome.fa
r1=$2/$1_R1_001.fastq.gz.cleaned.fastq.gz
#r2=$1_R2_001.fastq.gz                                                                                                                                                                        

mkdir ./$1
cd $1
#Alignment – Map to Reference
bwa mem -M -R '@RG\tID:sample_1\tLB:n12\tPL:ILLUMINA\tPM:HISEQ\tSM:n12' ${ref} ${r1} > $1_aligned_reads.sam

#Sort SAM file by coordinate, convert to BAM
java -jar ${PICARD_ROOT}/picard.jar SortSam INPUT=$1_aligned_reads.sam OUTPUT=$1_sorted_reads.bam SORT_ORDER=coordinate

#Collect Alignment & Insert Size Metrics
java -jar ${PICARD_ROOT}/picard.jar CollectAlignmentSummaryMetrics R=${ref} I=$1_sorted_reads.bam O=$1_alignment_metrics.txt
java -jar ${PICARD_ROOT}/picard.jar CollectInsertSizeMetrics INPUT=$1_sorted_reads.bam OUTPUT=$1_insert_metrics.txt HISTOGRAM_FILE=$1_insert_size_histogram.pdf
samtools depth -a $1_sorted_reads.bam > $1_depth_out.txt

#Mark Duplicates
java -jar ${PICARD_ROOT}/picard.jar MarkDuplicates INPUT=$1_sorted_reads.bam OUTPUT=$1_dedup_reads.bam METRICS_FILE=$1_metrics.txt

#Build BAM Index
java -jar ${PICARD_ROOT}/picard.jar BuildBamIndex INPUT=$1_dedup_reads.bam

#Create Realignment Targets
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref} -I $1_dedup_reads.bam -o $1_realignment_targets.list

#Realign Indels
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref} -I $1_dedup_reads.bam -targetIntervals $1_realignment_targets.list -o $1_realigned_reads.bam

#Call Variants
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${ref} -I $1_realigned_reads.bam -o $1_raw_variants.vcf

#Extract SNPs & Indels
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V $1_raw_variants.vcf -selectType SNP -o $1_raw_snps.vcf
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V $1_raw_variants.vcf -selectType INDEL -o $1_raw_indels.vcf

#Filter SNPs
#Note: SNPs which are ‘filtered out’ at this step will remain in the filtered_snps.vcf file, however they will be marked as ‘basic_snp_filter’, while SNPs which passed the filter will be marked as ‘PASS’
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref} -V $1_raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o $1_filtered_snps.vcf

#Filter Indels
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref} -V $1_raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o $1_filtered_indels.vcf

#Base Quality Score Recalibration (BQSR) #1
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref} -I $1_realigned_reads.bam -knownSites $1_filtered_snps.vcf -knownSites $1_filtered_indels.vcf -o $1_recal_data.table

#Base Quality Score Recalibration (BQSR) #2
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref} -I $1_realigned_reads.bam -knownSites $1_filtered_snps.vcf -knownSites $1_filtered_indels.vcf -BQSR $1_recal_data.table -o $1_post_recal_data.table

#Analyze Covariates
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${ref} -before $1_recal_data.table -after $1_post_recal_data.table -plots $1_recalibration_plots.pdf

#Apply BQSR
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T PrintReads -R ${ref} -I $1_realigned_reads.bam -BQSR $1_recal_data.table -o $1_recal_reads.bam

#Call Variants
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${ref} -I $1_recal_reads.bam -o $1_raw_variants_recal.vcf

#Extract SNPs & Indels
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V $1_raw_variants_recal.vcf -selectType SNP -o $1_raw_snps_recal.vcf
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V $1_raw_variants_recal.vcf -selectType INDEL -o $1_raw_indels_recal.vcf

#Filter SNPs
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref} -V $1_raw_snps_recal.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o $1_filtered_snps_final.vcf

#Filter Indels
java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref} -V $1_raw_indels_recal.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o $1_filtered_indels_recal.vcf

#Annotate SNPs and Predict Effects



#Compute Coverage Statistics
bedtools genomecov -bga -ibam $1_recal_reads.bam > $1_genomecov.bedgraph