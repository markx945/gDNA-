#!/bin/bash


# 显示帮助信息
function show_help {
    echo "Usage: $0 -i input_path -o output_path -f ref_genome -d dbsnp"
    echo
    echo "-i input_path    The path where to look for clean fastq file."
    echo "-o output_path   The path where to save the bam result. Usually the upper path of shell_script"
    echo "-f ref_genome    The path where to save the reference geneme file. fasta format"
    echo "-d dbsnp         The path where to save the known snp site file."
    echo
    echo "Example: $0 -i /path/to/cleanfastq -o /path/to/upper_path_of_shell_script -f /path/to/reference_genome -d /path/to/dbsnp_file"
    echo "Example2: $0 -i ../clean_data -o ../ -f ../other_ref/hg38.fa -d ../other_ref/dbsnp.vcf"
    exit 1
}

# 初始化输入和输出变量
input_path=""
output_path=""
ref_genome=""
dbsnp=""

# 解析命令行选项和参数
while getopts "hi:o:f:d:" opt; do
    case "$opt" in
    h)
        show_help
        ;;
    i)
        input_path=$OPTARG
        ;;
    o)
        output_path=$OPTARG
        ;;
    f)
        ref_genome=$OPTARG
        ;;
    d)
        dbsnp=$OPTARG
        ;;
    *)
        show_help
        ;;
    esac
done

# 检查是否提供了输入路径
if [ -z "$input_path" ]; then
    echo "Error: No input path provided."
    show_help
fi

# 读取路径下的所有文件夹并保存到数组中
SAMPLES=()
while IFS=  read -r -d $'\0'; do
    dir_name=$(basename "$REPLY")
    SAMPLES+=("$dir_name")
done < <(find "$input_path" -mindepth 1 -maxdepth 1 -type d -print0)


# 继续假设配置变量已经设置
# REF="/home/cfff_r2636/data/reference/hg38/genome/hg38.fa"
# DBSNP="/home/cfff_r2636/data/reference/hg38/dbsnp/dbsnp_146.hg38.vcf.gz"

MAPPING_OUTPUT_PATH="${output_path}mapping"
REF=${ref_genome}
DBSNP=${dbsnp}
QUALIMAP_OUTDIR="${output_path}qualimap"


for SAMPLE in "${SAMPLES[@]}"; do

    # bwa
    mkdir -p "${MAPPING_OUTPUT_PATH}/${SAMPLE}"
    bwa mem -M -R "@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA" -t 16 -K 10000000 "$REF" "${input_path}${SAMPLE}/${SAMPLE}_1.trimmed.fq.gz" "${input_path}${SAMPLE}/${SAMPLE}_2.trimmed.fq.gz" | \
    samtools view -@ 16 -b - | \
    samtools sort -@ 16 -o "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.sorted.bam" 

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
        -O "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal_data.table"

    # apply_bqsr
    gatk ApplyBQSR \
        -R "$REF" \
        -I "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.dedup.bam" \
        --bqsr-recal-file "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal_data.table" \
        -O "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal.bam"


    # generate bed file
    # 首先使用samtools查找高覆盖度区域
    samtools depth -a  ${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal.bam | awk '$3 > 10 {print $1"\t"$2-1"\t"$2}' > "${MAPPING_OUTPUT_PATH}/${SAMPLE}/high_coverage_10X.bed"

    # 使用bedtools合并接近的区域
    bedtools merge -i "${MAPPING_OUTPUT_PATH}/${SAMPLE}/high_coverage_10X.bed" > "${MAPPING_OUTPUT_PATH}/${SAMPLE}/merged_10X_coverage.bed"

    awk '{print $0 "\t0\t0\t0"}' "${MAPPING_OUTPUT_PATH}/${SAMPLE}/merged_10X_coverage.bed" > "${MAPPING_OUTPUT_PATH}/${SAMPLE}/merged_10X_coverage_6_line.bed"

    # qualimap
    QUALIMAP_SAMPLE_OUTDIR="${QUALIMAP_OUTDIR}/${SAMPLE}"
    mkdir -p "$QUALIMAP_SAMPLE_OUTDIR"
    qualimap bamqc \
        -bam "${MAPPING_OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.recal.bam" \
        -outformat PDF:HTML \
        -gff ${MAPPING_OUTPUT_PATH}/${SAMPLE}/merged_10X_coverage_6_line.bed \
        -nt 16 \
        -nr 500 \
        -nw 1500 \
        --java-mem-size=64G \
        -outdir "$QUALIMAP_SAMPLE_OUTDIR" \

done

echo "Bam file preprocess finished."