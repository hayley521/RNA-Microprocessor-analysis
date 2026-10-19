#!/bin/bash
#SBATCH -J GRO
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 

mkdir bw
for sample in R0-1 R0-4-2_R0-2 R0-4-2_R4-2 R4-1 R8_R8-1 R8_R8-2 R12_R12-1 R12_R12-2
do
#fastqc rawdata/${sample}_1.fastq.gz rawdata/${sample}_2.fastq.gz -t 8 -o fastqFileQC
source /home/xiegang/miniconda3/bin/activate
conda activate DNA 
bamCoverage -b dedup/${sample}.bam -o bw/${sample}.bw --normalizeUsing RPKM -p 18 --binSize 20
done
