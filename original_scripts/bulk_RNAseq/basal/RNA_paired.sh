source /home/xinhuil/anaconda3/bin/activate


for sample in CRGI_G13-1 CRGI_G13-2 CRGI_I13-1 CRGI_I13-2 CRGI_R1-1 CRGI_R1-2 CRGI_WT-1 CRGI_WT-2
do
fastqc rawdata/${sample}_1.fastq.gz rawdata/${sample}_2.fastq.gz -t 8 -o fastqFileQC
#conda activate trim-galore
#trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired /home/xinhuil/lyw-GROseq/260226/X101SC25071983-Z01-J205/01.RawData/${sample}/${sample}_1.fq.gz /home/xinhuil/lyw-GROseq/260226/X101SC25071983-Z01-J205/01.RawData/${sample}/${sample}_2.fq.gz --gzip -o ./clean --basename ${sample}
#conda activate RNA


#STAR --runThreadN 20 \
#--runMode alignReads \
#--genomeDir /home/xinhuil/reference/mouse \
#--readFilesCommand zcat --readFilesIn clean/${sample}_val_1.fq.gz clean/${sample}_val_2.fq.gz \
#--outFilterMultimapNmax 3 \
#--outFilterMismatchNoverLmax 0.04 \
#--outFilterType BySJout \
#--outSAMattrRGline ID:${sample} SM:${sample} PU:Illumina \
#--outSAMattributes All \
#--outSAMtype BAM SortedByCoordinate \
#--outFileNamePrefix ./align/${sample}_

done

#featureCounts -a /home/xinhuil/reference/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf -g gene_id -o rawCounts_O.tsv -T 20 -t exon -M -O -p --fraction \
 align/CRGI_G13-1_Aligned.sortedByCoord.out.bam \
 align/CRGI_G13-2_Aligned.sortedByCoord.out.bam \
 align/CRGI_I13-1_Aligned.sortedByCoord.out.bam \
 align/CRGI_I13-2_Aligned.sortedByCoord.out.bam \
 align/CRGI_R1-1_Aligned.sortedByCoord.out.bam \
 align/CRGI_R1-2_Aligned.sortedByCoord.out.bam \
 align/CRGI_WT-1_Aligned.sortedByCoord.out.bam \
 align/CRGI_WT-2_Aligned.sortedByCoord.out.bam 