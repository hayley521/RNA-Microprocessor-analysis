#!/bin/bash
#SBATCH -J rRNA
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log
#bowtie2-build mus_rRNA.fasta mus_rRNA
for sample in G4_G4-1 G4_G4-2 G8-1 G8-2 G12-1 G12-2
do
bowtie2 -x /home/xiegang/data/WKW/mus_rRNA --un-gz /home/xiegang/data/WKW/250211/clean/${sample}_rmrRNA.fastq.gz -U /home/xiegang/data/WKW/250211/clean/${sample}_trimmed.fq.gz -p 18 -S ${sample}_rRNA.sam
done
