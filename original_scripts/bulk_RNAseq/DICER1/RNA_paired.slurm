#!/bin/bash
#SBATCH -J DLF
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 

source /home/xiegang/miniconda3/bin/activate
conda activate base
export PATH="/home/xiegang/apps/TrimGalore-0.6.7:$PATH"
export PATH="/home/xiegang/apps/gatk-4.2.6.1:$PATH"

for sample in I1_I13-0h-1 I1_I13-1h-1 I1_I13-2h-1 I1_I13-4h-1 I2_I13-8h-1 I2_I13-12h-1 I2_I13-24h-1 I2_I13-72h-1 I3_I13-0h-2 I3_I13-1h-2 I3_I13-2h-2 I3_I13-4h-2 I4_I13-8h-2 I4_I13-12h-2 I4_I13-24h-2 I4_I13-72h-2
do

#fastqc rawdata/${sample}_1.fastq.gz rawdata/${sample}_2.fastq.gz -t 8 -o fastqFileQC

trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired /home/xiegang/data/WKW/250531/rawdata/X101SC24106259-Z01-J015/01.RawData/${sample}/${sample}_1.fq.gz /home/xiegang/data/WKW/250531/rawdata/X101SC24106259-Z01-J015/01.RawData/${sample}/${sample}_2.fq.gz --gzip -o ./clean --basename ${sample}

STAR --runThreadN 20 \
--runMode alignReads \
--genomeDir /home/xiegang/data/Reference/GRCm38/index_vM25 \
--readFilesCommand zcat --readFilesIn clean/${sample}_val_1.fq.gz clean/${sample}_val_2.fq.gz  \
--outFilterMultimapNmax 3 \
--outFilterMismatchNoverLmax 0.04 \
--outFilterType BySJout \
--outSAMattrRGline ID:${sample} SM:${sample} PU:Illumina \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./align/${sample}_

done

featureCounts -a /home/xiegang/data/Reference/GRCm38/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf -g gene_id -o rawCounts_O.tsv -T 20 -t exon -M -O -p --fraction \
 /home/xiegang/data/WKW/250531/align/I1_I13-0h-1_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I1_I13-1h-1_Aligned.sortedByCoord.out.bam \
  /home/xiegang/data/WKW/250531/align/I1_I13-2h-1_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I1_I13-4h-1_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I2_I13-8h-1_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I2_I13-12h-1_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I2_I13-24h-1_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I2_I13-72h-1_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I3_I13-0h-2_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I3_I13-1h-2_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I3_I13-2h-2_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I3_I13-4h-2_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I4_I13-8h-2_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I4_I13-12h-2_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I4_I13-24h-2_Aligned.sortedByCoord.out.bam \
 /home/xiegang/data/WKW/250531/align/I4_I13-72h-2_Aligned.sortedByCoord.out.bam 
