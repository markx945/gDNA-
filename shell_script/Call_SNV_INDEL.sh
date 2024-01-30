#! /bin/bash

set -euo pipefail

# 显示帮助信息
function show_help {
    echo "Usage: $0 -i input_path [-o output_file]"
    echo
    echo "-i normal_bam_file    The path where to look for normal bam file."
    echo "-I tumor_bam_file    The path where to look for tumor bam file."
    echo "-o output_path   The path where to save the SNV and INDEL result. Usually the upper path of shell_script"
    echo
    echo "Example: $0 -i /path/to/normal_bam_file -I /path/to/tumor_bam_file -o /path/to/upper_path_of_shell_script"
    echo
    exit 1
}

# 初始化输入和输出变量
normal_bam_file=""
tumor_bam_file=""
output_path=""

# 解析命令行选项和参数
while getopts "hi:I:o:" opt; do
    case "$opt" in
    h)
        show_help
        ;;
    i)
        normal_bam_file=$OPTARG
        ;;
    I)
        tumor_bam_file=$OPTARG
        ;;
    o)
        output_path=$OPTARG
        ;;
    *)
        show_help
        ;;
    esac
done


ref_fasta="/home/cfff_r2636/data/reference/hg38/genome/hg38.fa"
# ref_bed="/home/cfff_r2636/data/gDNA_231114/merged_10X_coverage.bed.gz"
out_dir="${output_path}vcf"
nt="16"


normal_bam=${normal_bam_file}
tumor_bam=${tumor_bam_file}

Normal_full_filename=$(basename "$normal_bam_file")
Normal_sample=$(echo "$Normal_full_filename" | cut -d '.' -f 1)
Tumor_full_filename=$(basename "$tumor_bam_file")
Tumor_sample=$(echo "$Tumor_full_filename" | cut -d '.' -f 1)

mkdir -p ${out_dir}/${Tumor_sample}

# 运行GATK Mutect2
gatk Mutect2 \
    -R "${ref_fasta}" \
    -I "${tumor_bam}" \
    -I "${normal_bam}" \
    --native-pair-hmm-threads ${nt} \
    -normal ${Normal_sample} \
    --tumor-sample ${Tumor_sample} \
    -O ${out_dir}/${Tumor_sample}/${Tumor_sample}.vcf

## 过滤vcf文件
gatk FilterMutectCalls \
    -V ${out_dir}/${Tumor_sample}/${Tumor_sample}.vcf \
    -R ${ref_fasta} \
    -O ${out_dir}/${Tumor_sample}/${Tumor_sample}.filter.vcf


### 注释文件
annotation_out_dir="${output_path}annovar"

dir_to_annovar="/home/cfff_r2636/data/software/annovar/annovar"
dir_db="/home/cfff_r2636/data/software/annovar/annovar"

mkdir -p ${annotation_out_dir}/${Tumor_sample}

perl ${dir_to_annovar}/table_annovar.pl ${out_dir}/${Tumor_sample}/${Tumor_sample}.filter.vcf ${dir_db}/humandb -buildver hg38  -out ${annotation_out_dir}/${Tumor_sample}/${Tumor_sample}  -remove  -protocol refGene  -operation g  -nastring .  -vcfinput  -thread 4

