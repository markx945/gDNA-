#! /bin/bash

# SAMPLES=(AF-D-1_combined AF-D-2_combined AF-P-1_combined AF-P-2_combined M8_1_combined M8_2_combined AF-S-1_combined AF-S-2_combined)
# SAMPLES=(HCC1395_combined HCC13951_combined)
SAMPLES=(HCC1395BL_combined HCC1395BL1_combined)
MAPPING_OUTPUT_PATH="/home/cfff_r2636/data/gDNA_231213/mapping"
QUALIMAP_OUTDIR="/home/cfff_r2636/data/gDNA_231213/qualimap_output"
# bed_file="/home/cfff_r2636/data/gDNA_231128/SureSelectXT_HS_V8_S33266340_hg38_Regions.bed"
# bed_file="/home/cfff_r2636/data/gDNA_231128/WES.bed"

bed_file="/home/cfff_r2636/data/gDNA_231213/merged_10X_coverage_6_line.bed"

for SAMPLE in "${SAMPLES[@]}"; do

    QUALIMAP_SAMPLE_OUTDIR="${QUALIMAP_OUTDIR}/${SAMPLE}_bed"
    mkdir -p "$QUALIMAP_SAMPLE_OUTDIR"
    qualimap bamqc \
        -bam "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal.bam" \
        -outformat PDF:HTML \
        -gff ${bed_file} \
        --java-mem-size=16G \
        -nt 16 \
        -nr 500 \
        -nw 1500 \
        -outdir "$QUALIMAP_SAMPLE_OUTDIR" \

done