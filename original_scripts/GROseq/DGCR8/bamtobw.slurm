#!/bin/bash
#SBATCH -J GRO
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 
mkdir bw
for sample in G0_G0-1 G0_G0-2 G4_G4-1 G4_G4-2 G8-1 G8-2 G12-1 G12-2
do
#fastqc rawdata/${sample}_1.fastq.gz rawdata/${sample}_2.fastq.gz -t 8 -o fastqFileQC
source /home/xiegang/miniconda3/bin/activate
conda activate DNA
bamCoverage -b dedup/${sample}.bam -o bw/${sample}.bw --normalizeUsing RPKM -p 18 --binSize 20
done

