#!/bin/bash

# SAMPLES=(AF_1_combined AF_2_combined HCC1395_combined HCC1395BL_combined M8_1_combined M8_2_combined P5_1_combined P5_2_combined)
SAMPLES=(AF_1_combined AF_2_combined)

# 定义变量
GENOME_LIB_DIR="/home/cfff_r2636/data/test_gRNA_ref/SV/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir" # 更改为你的genome library目录的实际路径
# LEFT_FQ="/home/cfff_r2636/data/test_gRNA_ref/fastp/B34_AF_2/E200011455_L01_B34_AF_2_1.trimmed.fq.gz"              # 更改为你的left (forward) reads文件的实际路径
# RIGHT_FQ="/home/cfff_r2636/data/test_gRNA_ref/fastp/B34_AF_2/E200011455_L01_B34_AF_2_2.trimmed.fq.gz"             # 更改为你的right (reverse) reads文件的实际路径
OUTPUT_DIR="/home/cfff_r2636/data/gDNA_231114/SV"   # 更改为你想要输出结果的目录的实际路径



for SAMPLE in "${SAMPLES[@]}"; do
    star_fusion_dir="/home/cfff_r2636/data/software/STAR-Fusion/STAR-Fusion"
    # 检查输出目录是否存在，如果不存在，则创建
    mkdir -p "$OUTPUT_DIR"/${SAMPLE}

    LEFT_FQ="/home/cfff_r2636/data/gDNA_231114/fastp/${SAMPLE}/${SAMPLE}_1.trimmed.fq.gz"            
    RIGHT_FQ="/home/cfff_r2636/data/gDNA_231114/fastp/${SAMPLE}/${SAMPLE}_2.trimmed.fq.gz" 

    # 运行STAR-Fusion
    ${star_fusion_dir} \
        --genome_lib_dir "$GENOME_LIB_DIR" \
        --left_fq "$LEFT_FQ" \
        --right_fq "$RIGHT_FQ" \
        --output_dir "$OUTPUT_DIR"/${SAMPLE}\
        --CPU 16

done