# build reference sequences
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz
seqkit seq Mus_musculus.GRCm38.ncrna.fa -n | grep -v lncRNA  | grep -v lincRNA > filter_list
seqkit grep -n -f filter_list Mus_musculus.GRCm38.ncrna.fa -o Mus_musculus.GRCm38.small.ncRNA.fa
bowtie-build Mus_musculus.GRCm38.small.ncRNA.fa Mus_musculus.GRCm38.small.ncRNA.fa


for sample in fastq/*
do
sample_name=`basename $sample`
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 ${sample}/${sample_name}.fq.gz -o ./clean --basename ${sample_name:2} --dont_gzip
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/${sample_name:2}_trimmed.fq  -S align/${sample_name:2}.sam --un align/${sample_name:2}.unaligned.fq --allow-contain -p 40 
samtools view -h -F 4 align/${sample_name:2}.sam -@ 20 | samtools sort -@ 20 - -o align/${sample_name:2}.sorted.sam
samtools view align/${sample_name:2}.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > ${sample_name:2}.counts.txt
done

for sample in fastq/*
do
sample_name=`basename $sample`
pigz -d clean/${sample_name:2}_trimmed.fq.gz
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/mature_mmu_spike.fa \
 -q clean/${sample_name:2}_trimmed.fq  -S align/${sample_name:2}.sam --un align/${sample_name:2}.unaligned.fq --allow-contain -p 40 
samtools view -h -F 4 align/${sample_name:2}.sam -@ 20 | samtools sort -@ 20 - -o align/${sample_name:2}.sorted.sam
samtools view align/${sample_name:2}.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > ${sample_name:2}.maturecounts.txt
pigz clean/${sample_name:2}_trimmed.fq
done

# trim adaptor
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J037/01.RawData/G5mi_G5C/G5mi_G5C.fq.gz -o ./clean --basename G5C --dont_gzip
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J037/01.RawData/G5mi_G5I/G5mi_G5I.fq.gz -o ./clean --basename G5I --dont_gzip
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J063/01.RawData/G5mi-2_G5C-2/G5mi-2_G5C-2.fq.gz -o ./clean --basename G5C-2 --dont_gzip
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J063/01.RawData/G5mi-2_G5I-2/G5mi-2_G5I-2.fq.gz -o ./clean --basename G5I-2 --dont_gzip

# map to small ncRNA
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/G5C_trimmed.fq  -S align/G5C.sam --un align/G5C.unaligned.fq --allow-contain -p 40 
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/G5I_trimmed.fq  -S align/G5I.sam --un align/G5I.unaligned.fq --allow-contain -p 40 
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/G5C-2_trimmed.fq  -S align/G5C-2.sam --un align/G5C-2.unaligned.fq --allow-contain -p 40 
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/G5I-2_trimmed.fq  -S align/G5I-2.sam --un align/G5I-2.unaligned.fq --allow-contain -p 40 

# map to genome
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/G5C.unaligned.fq -S align/G5C_genome.sam --allow-contain -p 40 
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/G5I.unaligned.fq -S align/G5I_genome.sam --allow-contain -p 40 
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/G5C-2.unaligned.fq -S align/G5C-2_genome.sam --allow-contain -p 40 
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/G5I-2.unaligned.fq -S align/G5I-2_genome.sam --allow-contain -p 40 

# sort and count
samtools view -h -F 4 align/G5C.sam -@ 20 | samtools sort -@ 20 - -o align/G5C.sorted.sam
samtools view align/G5C.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > G5C.counts.txt

samtools view -h -F 4 align/G5I.sam -@ 20 | samtools sort -@ 20 - -o align/G5I.sorted.sam
samtools view align/G5I.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > G5I.counts.txt

samtools view -h -F 4 align/G5C-2.sam -@ 20 | samtools sort -@ 20 - -o align/G5C-2.sorted.sam
samtools view align/G5C-2.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > G5C-2.counts.txt

samtools view -h -F 4 align/G5I-2.sam -@ 20 | samtools sort -@ 20 - -o align/G5I-2.sorted.sam
samtools view align/G5I-2.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > G5I-2.counts.txt

# trim adaptor
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J033/01.RawData/R0/R0.fq.gz -o ./clean --basename R0 --dont_gzip
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J033/01.RawData/R4/R4.fq.gz -o ./clean --basename R4 --dont_gzip
 trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J033/01.RawData/R8/R8.fq.gz -o ./clean --basename R8 --dont_gzip
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
 X102SC23111973-Z01-J033/01.RawData/R12/R12.fq.gz -o ./clean --basename R12 --dont_gzip

# map to small ncRNA
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/R0_trimmed.fq  -S align/R0.sam --un align/R0.unaligned.fq --allow-contain -p 40 
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/R4_trimmed.fq  -S align/R4.sam --un align/R4.unaligned.fq --allow-contain -p 40 
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/R8_trimmed.fq  -S align/R8.sam --un align/R8.unaligned.fq --allow-contain -p 40 
bowtie -a --best --strata -n 1  -x /wanglab1/gangx/data/LYW/ref/Mus_musculus.GRCm38.small.ncRNA.fa \
 -q clean/R12_trimmed.fq  -S align/R12.sam --un align/R12.unaligned.fq --allow-contain -p 40 

# map to genome
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/R0.unaligned.fq -S align/R0_genome.sam --allow-contain -p 40 
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/R4.unaligned.fq -S align/R4_genome.sam --allow-contain -p 40 
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/R8.unaligned.fq -S align/R8_genome.sam --allow-contain -p 40 
bowtie -a --best --strata -n 1 -m 10 -x /wanglab1/gangx/data/Reference/10090/GRCm38.p6/GRCm38.p6.genome.fa \
 -q align/R12.unaligned.fq -S align/R12_genome.sam --allow-contain -p 40 

# sort and count
samtools view -h -F 4 align/R0.sam -@ 20 | samtools sort -@ 20 - -o align/R0.sorted.sam
samtools view align/R0.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > R0.counts.txt

samtools view -h -F 4 align/R4.sam -@ 20 | samtools sort -@ 20 - -o align/R4.sorted.sam
samtools view align/R4.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > R4.counts.txt

samtools view -h -F 4 align/R8.sam -@ 20 | samtools sort -@ 20 - -o align/R8.sorted.sam
samtools view align/R8.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > R8.counts.txt

samtools view -h -F 4 align/R12.sam -@ 20 | samtools sort -@ 20 - -o align/R12.sorted.sam
samtools view align/R12.sorted.sam | awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]}END{for(i in x) print i, x[i]}' | sort  -gr -k 2 > R12.counts.txt