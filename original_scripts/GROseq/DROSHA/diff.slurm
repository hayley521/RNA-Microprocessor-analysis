#!/bin/bash
#SBATCH -J GRO
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 
analyzeRepeats.pl rna mm10 -count genes -condenseGenes -raw -d /home/xiegang/data/WKW/250224/HomerTag/R0-1 /home/xiegang/data/WKW/250224/HomerTag/R0-4-2_R0-2 /home/xiegang/data/WKW/250224/HomerTag/R0-4-2_R4-2 /home/xiegang/data/WKW/250224/HomerTag/R4-1 /home/xiegang/data/WKW/250224/HomerTag/R8_R8-1 /home/xiegang/data/WKW/250224/HomerTag/R8_R8-2 /home/xiegang/data/WKW/250224/HomerTag/R12_R12-1 /home/xiegang/data/WKW/250224/HomerTag/R12_R12-2  > countTable.txt
getDiffExpression.pl countTable.txt R0 R0 R4 R4 R8 R8 R12 R12 > diffOutput.txt