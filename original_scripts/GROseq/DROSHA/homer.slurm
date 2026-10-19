#!/bin/bash
#SBATCH -J GRO
#SBATCH -p c20m128 
#SBATCH -n 10
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 
#mkdir HomerTag
for sample in R0-1 R0-4-2_R0-2 R0-4-2_R4-2 R4-1 R8_R8-1 R8_R8-2 R12_R12-1 R12_R12-2
do
makeTagDirectory HomerTag/${sample} align/${sample}_Aligned.sorted.bam -genome mm10
done

#export PATH=/home/xiegang/data/WKW/250224/bin:$PATH
analyzeRepeats.pl rna mm10 -count genes -condenseGenes  -norm 1e6 -d HomerTag/R0-1 HomerTag/R0-4-2_R0-2 HomerTag/R4-1 HomerTag/R0-4-2_R4-2 HomerTag/R8_R8-1 HomerTag/R8_R8-2 HomerTag/R12_R12-1 HomerTag/R12_R12-2 > homer_genes_norm.txt
analyzeRepeats.pl rna mm10 -count pausing -condenseGenes -norm 1e6 -d HomerTag/R0-1 HomerTag/R0-4-2_R0-2 HomerTag/R4-1 HomerTag/R0-4-2_R4-2 HomerTag/R8_R8-1 HomerTag/R8_R8-2 HomerTag/R12_R12-1 HomerTag/R12_R12-2 > homer_pausing_norm.txt 
analyzeRepeats.pl rna mm10 -count genes -condenseGenes -raw -d /home/xiegang/data/WKW/250224/HomerTag/R0-1 /home/xiegang/data/WKW/250224/HomerTag/R0-4-2_R0-2 /home/xiegang/data/WKW/250224/HomerTag/R0-4-2_R4-2 /home/xiegang/data/WKW/250224/HomerTag/R4-1 /home/xiegang/data/WKW/250224/HomerTag/R8_R8-1 /home/xiegang/data/WKW/250224/HomerTag/R8_R8-2 /home/xiegang/data/WKW/250224/HomerTag/R12_R12-1 /home/xiegang/data/WKW/250224/HomerTag/R12_R12-2  > countTable.txt
