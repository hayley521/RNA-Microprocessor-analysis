#!/bin/bash
#SBATCH -J GRO
#SBATCH -p c20m128 
#SBATCH -n 10
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 

for sample in G0_G0-1 G0_G0-2 G4_G4-1 G4_G4-2 G8-1 G8-2 G12-1 G12-2
do
#fastqc rawdata/${sample}_1.fastq.gz rawdata/${sample}_2.fastq.gz -t 8 -o fastqFileQC
source /home/xiegang/miniconda3/bin/activate
conda activate base
export PATH="/home/xiegang/apps/TrimGalore-0.6.7:$PATH"
export PATH="/home/xiegang/apps/gatk-4.2.6.1:$PATH"
#trim_galore -j 10 -q 20 --clip_R1 3 --phred33 --stringency 3 --length 20 -e 0.1 /home/xiegang/data/WKW/250211/X101SC24090344-Z01-J082/01.RawData/${sample}/${sample}.fq.gz --gzip -o ./clean --basename ${sample}
#trim_galore --polyA clean/${sample}_trimmed.fq.gz --gzip
STAR --runThreadN 8 \
--runMode alignReads \
--genomeDir /data1/xiegang/data/Reference/GRCm38/index_vM25/ \
--readFilesCommand zcat --readFilesIn ${sample}_trimmed_trimmed.fq.gz  \
--outFilterMultimapNmax 3 \
--outFilterMismatchNoverLmax 0.04 \
--outFilterType BySJout \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./align/${sample}_
samtools sort -@ 10 align/${sample}_Aligned.sortedByCoord.out.bam -o ./align/${sample}_Aligned.sorted.bam 
samtools index -@ 10 ./align/${sample}_Aligned.sorted.bam
picard MarkDuplicates -I align/${sample}_Aligned.sorted.bam -O dedup/${sample}.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M dedup/${sample}_MarkDuplicates_output.metrics 

conda activate DNA
bamCoverage -b dedup/${sample}.bam -o align/${sample}_forward.bw --normalizeUsing RPKM -p 8 --binSize 20 --filterRNAstrand forward
bamCoverage -b dedup/${sample}.bam -o align/${sample}_reverse.bw --normalizeUsing RPKM -p 8 --binSize 20 --filterRNAstrand reverse 
conda activate macs2
macs2 callpeak -t dedup/${sample}.bam --keep-dup all --nolambda --nomodel -B --extsize 150  -g hs --outdir peak/ -n ${sample}

done

