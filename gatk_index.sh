#!/bin/bash                                          

#SBATCH --job-name=index
##SBATCH --nodes=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
##SBATCH --gres=gpu:1
##SBATCH --partition=gpu4_medium
#SBATCH --partition=cpu_long
#SBATCH --error=/gpfs/data/proteomics/projects/Sunny/jingchuan/yeast_genome/err_out/%x_%j.err
#SBATCH --output=/gpfs/data/proteomics/projects/Sunny/jingchuan/yeast_genome/err_out/%x_%j.out
##SBATCH --dependency=afterany:job_id

module load bwa/0.7.17
module load picard/2.18.11
module load samtools/1.3
module load gatk/3.6.0
module load bowtie2

bowtie2-build $1.fa $1

bwa index $1.fa

java -jar ${PICARD_ROOT}/libs/picard.jar CreateSequenceDictionary R=$1.fa O=$1.dict

samtools faidx $1.fa
