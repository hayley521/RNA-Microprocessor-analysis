#!/bin/bash
#SBATCH -J DLF
#SBATCH -p c20m128 
#SBATCH -n 20
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log 
# 设置路径
GTF=/home/xiegang/data/Reference/GRCm38/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf
BAM_DIR=/home/xiegang/data/WKW/250531/align
OUTPUT_DIR=rmats_results
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
mkdir -p $OUTPUT_DIR

# 定义时间点（排除0h）
for time in 1h 2h 4h 8h 12h 24h 72h; do
    echo "Processing ${time} vs 0h..."
    
    # 创建该时间点的bam列表
    # 注意：0h有两个重复（I1和I3），每个时间点也有两个重复
    cat > ${OUTPUT_DIR}/bams_0h.txt << EOF
${BAM_DIR}/I1_I13-0h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I3_I13-0h-2_Aligned.sortedByCoord.out.bam
EOF
    
    # 根据时间点选择对应的bam文件
    case $time in
        1h)
            cat > ${OUTPUT_DIR}/bams_${time}.txt << EOF
${BAM_DIR}/I1_I13-1h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I3_I13-1h-2_Aligned.sortedByCoord.out.bam
EOF
            ;;
        2h)
            cat > ${OUTPUT_DIR}/bams_${time}.txt << EOF
${BAM_DIR}/I1_I13-2h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I3_I13-2h-2_Aligned.sortedByCoord.out.bam
EOF
            ;;
        4h)
            cat > ${OUTPUT_DIR}/bams_${time}.txt << EOF
${BAM_DIR}/I1_I13-4h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I3_I13-4h-2_Aligned.sortedByCoord.out.bam
EOF
            ;;
        8h)
            cat > ${OUTPUT_DIR}/bams_${time}.txt << EOF
${BAM_DIR}/I2_I13-8h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I4_I13-8h-2_Aligned.sortedByCoord.out.bam
EOF
            ;;
        12h)
            cat > ${OUTPUT_DIR}/bams_${time}.txt << EOF
${BAM_DIR}/I2_I13-12h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I4_I13-12h-2_Aligned.sortedByCoord.out.bam
EOF
            ;;
        24h)
            cat > ${OUTPUT_DIR}/bams_${time}.txt << EOF
${BAM_DIR}/I2_I13-24h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I4_I13-24h-2_Aligned.sortedByCoord.out.bam
EOF
            ;;
        72h)
            cat > ${OUTPUT_DIR}/bams_${time}.txt << EOF
${BAM_DIR}/I2_I13-72h-1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/I4_I13-72h-2_Aligned.sortedByCoord.out.bam
EOF
            ;;
    esac
    
    # 运行rMATS
    rmats.py \
        --b1 ${OUTPUT_DIR}/bams_0h.txt \
        --b2 ${OUTPUT_DIR}/bams_${time}.txt \
        --gtf $GTF \
        --od ${OUTPUT_DIR}/${time}_vs_0h \
        --tmp ${OUTPUT_DIR}/tmp_${time} \
    -t paired \
    --readLength 150 \
    --nthread 18 \
    --libType fr-unstranded \
    --cstat 0.0001 \
    --novelSS
    
    echo "Finished ${time} vs 0h"
done

echo "All timepoints completed!"