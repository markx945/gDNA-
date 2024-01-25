#!/bin/bash

source /usr/conda/etc/profile.d/conda.sh

# fastp,mapping
conda activate mapping

# SAMPLES=(HCC1395_combined HCC1395BL_combined HCC13951_combined HCC1395BL1_combined)
SAMPLES=(HCC1395_combined HCC1395BL_combined)

REF="/home/cfff_r2636/data/reference/hg38/genome/hg38.fa"

TRIMMED_OUTPUT_PATH="/home/cfff_r2636/data/gDNA_231213/fastp"
MAPPING_OUTPUT_PATH="/home/cfff_r2636/data/gDNA_231213/mapping"


for SAMPLE in "${SAMPLES[@]}"; do
    # 根据SAMPLE变量设置相应的输入文件路径
    FASTQ_1="/home/cfff_r2636/data/gDNA_231213/FASTQ1/${SAMPLE}_R1.fastq.gz"
    FASTQ_2="/home/cfff_r2636/data/gDNA_231213/FASTQ1/${SAMPLE}_R2.fastq.gz"
    

    # fastp
    mkdir -p "${TRIMMED_OUTPUT_PATH}/${SAMPLE}"
    fastp --thread 16 -i ${FASTQ_1} -I ${FASTQ_2} -o "${TRIMMED_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}_1.trimmed.fq.gz" -O "${TRIMMED_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}_2.trimmed.fq.gz" -h ${TRIMMED_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.html

    # bwa
    mkdir -p "${MAPPING_OUTPUT_PATH}/${SAMPLE}"
    bwa mem -M -R "@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA" -t 16 -K 10000000 "$REF" "${TRIMMED_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}_1.trimmed.fq.gz" "${TRIMMED_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}_2.trimmed.fq.gz" | \
    samtools view -@ 16 -b - | \
    samtools sort -@ 16 -o "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.sorted.bam"

done

# 继续假设配置变量已经设置
BQSR_OUTPUT_PATH="/home/cfff_r2636/data/gDNA_231213/mapping"
# REF="/home/cfff_r2636/data/reference/hg38/genome/hg38.fa"
DBSNP="/home/cfff_r2636/data/reference/hg38/dbsnp/dbsnp_146.hg38.vcf.gz"
# GFF="/path/to/annotation.gff"
QUALIMAP_OUTDIR="/home/cfff_r2636/data/gDNA_231213/qualimap_output"

conda activate variant_calling

for SAMPLE in "${SAMPLES[@]}"; do

    picard MarkDuplicates \
        I="${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.sorted.bam" \
        O="${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.dedup.bam" \
        M="${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.metrics.txt" \
        REMOVE_DUPLICATES=true
    samtools index "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.dedup.bam"

    # base_recalibrator
    gatk BaseRecalibrator \
        -I "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.dedup.bam" \
        -R "$REF" \
        --known-sites "$DBSNP" \
        -O "${BQSR_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal_data.table"

    # apply_bqsr
    gatk ApplyBQSR \
        -R "$REF" \
        -I "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.dedup.bam" \
        --bqsr-recal-file "${BQSR_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal_data.table" \
        -O "${BQSR_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal.bam"

    # # qualimap
    # QUALIMAP_SAMPLE_OUTDIR="${QUALIMAP_OUTDIR}/${SAMPLE}_qualimap"
    # mkdir -p "$QUALIMAP_SAMPLE_OUTDIR"
    # qualimap bamqc \
    #     -bam "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal.bam" \
    #     -outformat PDF:HTML \
    #     -nt 16 \
    #     -nr 500 \
    #     -nw 1500 \
    #     --java-mem-size=64G \
    #     -outdir "$QUALIMAP_SAMPLE_OUTDIR" \

done
