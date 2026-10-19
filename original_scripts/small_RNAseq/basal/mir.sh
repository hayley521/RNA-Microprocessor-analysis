source /home/xinhuil/anaconda3/bin/activate
conda activate RNA
bowtie-build /home/xinhuil/reference/microRNA/mmu_mature_dna.fa /home/xinhuil/reference/microRNA/miRBase_mmu_bowtie -p 20
for sample in G_G13-1 G_G13-2 I_I13-1 I_I13-2 R_R1-1 R_R1-2 WT_WT-1 WT_WT-2
do
conda activate trim-galore
# 使用Trim Galore进行修剪
#trim_galore --small_rna --length 15 /home/xinhuil/lyw-GROseq/260227/X101SC25071983-Z01-J212/01.RawData/${sample}/${sample}.fq.gz \
#          -o clean/ --basename MiR-${sample} -j 20 --dont_gzip
conda activate RNA
bowtie -n 1 -a --best --strata /home/xinhuil/reference/microRNA/miRBase_mmu_bowtie \
       clean/MiR-${sample}_trimmed.fq \
       -S align/MiR-${sample}.sam -p 20

# SAM转BAM并排序
samtools view -bS align/MiR-${sample}.sam | samtools sort -o MiR-${sample}_sorted.bam

done
for bam in *.bam; do
    samtools index "$bam"
    samtools idxstats "$bam" > "${bam%.bam}.idxstats.txt"
done

# 提取第一个文件的前两列（miRNA和length）
cut -f1,2 MiR-WT_WT-1_sorted.idxstats.txt > cols1_2.txt

# 循环提取每个文件的第三列并保存为临时文件
for i in WT_WT-1 WT_WT-2 G_G13-1 G_G13-2 I_I13-1 I_I13-2 R_R1-1 R_R1-2; do
    cut -f3 MiR-${i}_sorted.idxstats.txt > temp_col${i}.txt
done

# 合并所有临时文件
paste cols1_2.txt temp_colWT_WT-1.txt temp_colWT_WT-2.txt temp_colG_G13-1.txt temp_colG_G13-2.txt temp_colI_I13-1.txt temp_colI_I13-2.txt temp_colR_R1-1.txt temp_colR_R1-2.txt > merged_miRNA_mature_2.txt
# 清理临时文件
rm cols1_2.txt temp_col*.txt
