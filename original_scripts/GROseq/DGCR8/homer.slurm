#!/bin/bash
#SBATCH -J GRO
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 

for sample in G0_G0-1 G0_G0-2 G4_G4-1 G4_G4-2 G8-1 G8-2 G12-1 G12-2
do
makeTagDirectory HomerTag/${sample} align/${sample}_Aligned.sorted.bam -genome mm10
done

analyzeRepeats.pl rna mm10 -count genes -raw -condenseGenes -d HomerTag/G0_G0-1 HomerTag/G0_G0-2 HomerTag/G4_G4-1 HomerTag/G4_G4-2 HomerTag/G8-1 HomerTag/G8-2 HomerTag/G12-1 HomerTag/G12-2 > countTable.txt
analyzeRepeats.pl rna mm10 -count pausing -condenseGenes -norm 1e6 -d HomerTag/G0_G0-1 HomerTag/G0_G0-2 HomerTag/G4_G4-1 HomerTag/G4_G4-2 HomerTag/G8-1 HomerTag/G8-2 HomerTag/G12-1 HomerTag/G12-2 > homer_pausing_norm.txt 
analyzeRepeats.pl rna mm10 -count genes -condenseGenes -norm 1e6 -d HomerTag/G0_G0-1 HomerTag/G0_G0-2 HomerTag/G4_G4-1 HomerTag/G4_G4-2 HomerTag/G8-1 HomerTag/G8-2 HomerTag/G12-1 HomerTag/G12-2 > homer_genes_norm.txt 

