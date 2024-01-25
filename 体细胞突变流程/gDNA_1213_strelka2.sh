#! /bin/bash

source /usr/conda/etc/profile.d/conda.sh

conda activate python2

set -euo pipefail

# Normal_sample=(M8-1_combined M8-2_combined M8-1_combined M8-2_combined)
# Tumor_sample=(AF-D-1_combined AF-D-2_combined AF-P-1_combined AF-P-2_combined)

# Normal_sample=(M8-2_combined M8-1_combined M8-2_combined M8-1_combined M8-2_combined)
# Tumor_sample=(AF-D-2_combined AF-P-1_combined AF-P-2_combined AF-S-1_combined AF-S-2_combined)

Normal_sample=(HCC1395BL_combined HCC1395BL_combined)
Tumor_sample=(HCC1395_combined HCC13951_combined) 

sample_length=${#Normal_sample[@]}


dit_to_strelka="/home/cfff_r2636/data/software/strelka-2.9.2/strelka-2.9.2.centos6_x86_64/bin"
ref_fasta="/home/cfff_r2636/data/reference/hg38/genome/hg38.fa"
ref_bed="/home/cfff_r2636/data/gDNA_231128/merged_10X_coverage.bed.gz"
out_dir="/home/cfff_r2636/data/gDNA_231213/variant_calling"
nt="16"

set -exo pipefail

for ((i=0; i<$sample_length; i++)); do

    normal_bam="/home/cfff_r2636/data/gDNA_231213/mapping/"${Normal_sample[$i]}/${Normal_sample[$i]}.recal.bam
    tumor_bam="/home/cfff_r2636/data/gDNA_231213/mapping/"${Tumor_sample[$i]}/${Tumor_sample[$i]}.recal.bam

    # mkdir -p ${out_dir}/${Tumor_sample[$i]}
    mkdir -p ${out_dir}/${Tumor_sample[$i]}

    ${dit_to_strelka}/configureStrelkaSomaticWorkflow.py \
        --normalBam ${normal_bam} \
        --tumorBam ${tumor_bam} \
        --callRegions ${ref_bed} \
        --referenceFasta ${ref_fasta} \
        --runDir ${out_dir}/${Tumor_sample[$i]}

    # without bed file
    # ${dit_to_strelka}/configureStrelkaSomaticWorkflow.py \
    #     --normalBam ${normal_bam} \
    #     --tumorBam ${tumor_bam} \
    #     --referenceFasta ${ref_fasta} \
    #     --runDir ${out_dir}/${Tumor_sample[$i]}

    python2 ${out_dir}/${Tumor_sample[$i]}/runWorkflow.py -m local -j ${nt}
done


### filter_vcf
# dir_base=${out_dir}

# snvs_vcf="somatic.snvs.vcf.gz"
# indels_vcf="somatic.indels.vcf.gz"

# snvs_out_vcf="somatic.snvs.filter.vcf"
# indels_out_vcf="somatic.indels.filter.vcf"

# for dir in */ ; do
#     if [[ -d "$dir" ]]; then
#         folder_name=$(basename "$dir")
#         path_to_vcf=${dir_base}/${folder_name}/results/variants/

#         # zcat $path_to_vcf/$snvs_vcf | awk '!/^##/ {if (/^#/) {print $0} else if ($7 == "PASS") {print $0}}' > $path_to_vcf/$snvs_out_vcf
#         zcat $path_to_vcf/$snvs_vcf | awk '{if (/^#/) {print $0} else if ($7 == "PASS") {print $0}}' > $path_to_vcf/$snvs_out_vcf

#         # zcat $path_to_vcf/$indels_vcf | awk '!/^##/ {if (/^#/) {print $0} else if ($7 == "PASS") {print $0}}' > $path_to_vcf/$indels_out_vcf
#         zcat $path_to_vcf/$indels_vcf | awk '{if (/^#/) {print $0} else if ($7 == "PASS") {print $0}}' > $path_to_vcf/$indels_out_vcf
#         ## 合并snv和indel vcf文件
#         bcftools concat -a $path_to_vcf/$indels_vcf $path_to_vcf/$snvs_out_vcf -o ${path_to_vcf}/${folder_name}.strelka.pass.vcf
#         bcftools norm -m -both ${path_to_vcf}/${folder_name}.strelka.pass.vcf -o ${path_to_vcf}/${folder_name}.strelka.norm.pass.vcf
#         bcftools index -t ${path_to_vcf}/${folder_name}.strelka.norm.pass.vcf
        
#         ##准备annotation输入

#         awk '{ if ($1 ~ /^#/) print; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 }' ${path_to_vcf}/${folder_name}.strelka.norm.pass.vcf | sed '/^#/!s/ /\t/g' > ${path_to_vcf}/${folder_name}.strelka.annovar.vcf
#     fi
# done




### annotation

