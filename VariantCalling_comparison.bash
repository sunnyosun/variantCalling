#!/bin/bash

#SBATCH --job-name=VC_compare
##SBATCH --nodes=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
##SBATCH --gres=gpu:1
##SBATCH --partition=gpu4_medium
#SBATCH --partition=cpu_long
#SBATCH --error=/gpfs/data/proteomics/projects/Sunny/jingchuan/2018-07-24_syni/rnaseq/err_out/%x_%j.err
#SBATCH --output=/gpfs/data/proteomics/projects/Sunny/jingchuan/2018-07-24_syni/rnaseq/err_out/%x_%j.out
##SBATCH --dependency=afterany:job_id

module load r/3.4.3

# all the inputs need to include paths
# 1. VCF file 1
# 2. VCF file 2 (ctr)
# 3. reference file (gff)
# 4. output name (optional, default is out.txt)

ref=/gpfs/data/proteomics/projects/Sunny/jingchuan/yeast_genome/BY4741_genome.gff
vcf1_snp=$1_filtered_snps_final.vcf
vcf2_snp=$2_filtered_snps_final.vcf
vcf1_indel=$1_filtered_indels_recal.vcf
vcf2_indel=$2_filtered_indels_recal.vcf
outname=$1_vs_$2

Rscript /gpfs/data/proteomics/projects/Sunny/R/variant_calling_indels_forbash.R ${vcf1_indel} ${vcf2_indel} ${ref} ${outname}
Rscript /gpfs/data/proteomics/projects/Sunny/R/variant_calling_snps_forbash.R ${vcf1_snp} ${vcf2_snp} ${ref} ${outname}

