#!/bin/bash
# rMATS 分析: 对照 vs R12

/home/xinhuil/anaconda3/envs/rmats/bin/rmats.py \
    --b1 b1.txt \
    --b2 b2_R12.txt \
    --gtf /home/xinhuil/reference/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf \
    --od rmats_output_G12 \
    --tmp rmats_output_G12/tmp \
    -t paired \
    --readLength 150 \
    --nthread 8 \
    --libType fr-unstranded \
    --cstat 0.0001 \
    --novelSS

echo "分析完成: 对照 vs R12"
echo "结果保存在: rmats_output_G12"
