#!/bin/bash
#SBATCH -J GRO
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 


for sample in R0-1 R0-4-2_R0-2 R0-4-2_R4-2 R4-1 R8_R8-1 R8_R8-2 R12_R12-1 R12_R12-2
do
#fastqc rawdata/${sample}_1.fastq.gz rawdata/${sample}_2.fastq.gz -t 8 -o fastqFileQC
source /home/xiegang/miniconda3/bin/activate
conda activate base
#export PATH="/home/xiegang/apps/TrimGalore-0.6.7:$PATH"
export PATH="/home/xiegang/apps/gatk-4.2.6.1:$PATH"
trim_galore -j 20 -q 20 --clip_R1 3 --phred33 --stringency 3 --length 20 -e 0.1 /home/xiegang/data/WKW/250224/rawdata/X101SC24090344-Z01-J106/01.RawData/${sample}/${sample}.fq.gz --gzip -o clean --basename ${sample}
trim_galore --polyA clean/${sample}_trimmed.fq.gz --gzip
bowtie2 -x /home/xiegang/data/WKW/mus_rRNA --un-gz clean/${sample}_rmrRNA.fastq.gz -U ${sample}_trimmed.fq.gz -p 18 -S ${sample}_rRNA.sam
STAR --runThreadN 18 \
--runMode alignReads \
--genomeDir /data1/xiegang/data/Reference/GRCm38/index_vM25/ \
--readFilesCommand zcat \
--readFilesIn clean/${sample}_rmrRNA.fastq.gz \
--outFilterMultimapNmax 3 \
--outFilterMismatchNoverLmax 0.04 \
--outFilterType BySJout \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./align/${sample}_
samtools sort -@ 20 align/${sample}_Aligned.sortedByCoord.out.bam -o ./align/${sample}_Aligned.sorted.bam 
samtools index -@ 20 ./align/${sample}_Aligned.sorted.bam
picard MarkDuplicates -I align/${sample}_Aligned.sorted.bam -O dedup/${sample}.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M dedup/${sample}_MarkDuplicates_output.metrics 

conda activate DNA 
bamCoverage -b dedup/${sample}.bam -o align/${sample}_forward.bw --normalizeUsing RPKM -p 18 --binSize 20 --filterRNAstrand forward
bamCoverage -b dedup/${sample}.bam -o align/${sample}_reverse.bw --normalizeUsing RPKM -p 18 --binSize 20 --filterRNAstrand reverse 
conda activate macs2
macs2 callpeak -t dedup/${sample}.bam --keep-dup all --nolambda --nomodel -B --extsize 150  -g hs --outdir peak/ -n ${sample}

done


