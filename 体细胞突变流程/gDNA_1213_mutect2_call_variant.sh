#! /bin/bash

set -euo pipefail

# Normal_sample=(M8-1_combined M8-2_combined M8-1_combined M8-2_combined)
# Tumor_sample=(AF-D-1_combined AF-D-2_combined AF-P-1_combined AF-P-2_combined)


# Normal_sample=(M8-1_combined M8-2_combined)
# Tumor_sample=(AF-S-1_combined AF-S-2_combined)

# Normal_sample=(HCC1395BL_combined HCC1395BL_combined)
# Tumor_sample=(HCC1395_combined HCC13951_combined) 
Normal_sample=(HCC1395BL1_combined)
Tumor_sample=(HCC13951_combined) 

sample_length=${#Normal_sample[@]}


ref_fasta="/home/cfff_r2636/data/reference/hg38/genome/hg38.fa"
# ref_bed="/home/cfff_r2636/data/gDNA_231114/merged_10X_coverage.bed.gz"
out_dir="/home/cfff_r2636/data/gDNA_231213/variant_calling/mutect2"

nt="16"

set -exo pipefail

for ((i=0; i<$sample_length; i++)); do

    normal_bam="/home/cfff_r2636/data/gDNA_231213/mapping/"${Normal_sample[$i]}/${Normal_sample[$i]}.recal.bam
    tumor_bam="/home/cfff_r2636/data/gDNA_231213/mapping/"${Tumor_sample[$i]}/${Tumor_sample[$i]}.recal.bam


    # mkdir -p ${out_dir}/${Tumor_sample[$i]}
    mkdir -p ${out_dir}/${Tumor_sample[$i]}

    # 运行GATK Mutect2
    gatk Mutect2 \
    -R "${ref_fasta}" \
    -I "${tumor_bam}" \
    -I "${normal_bam}" \
    --native-pair-hmm-threads ${nt} \
    -normal ${Normal_sample[$i]} \
    --tumor-sample ${Tumor_sample[$i]} \
    -O ${out_dir}/${Tumor_sample[$i]}/${Tumor_sample[$i]}.vcf

    # gatk VariantFiltration \
    #                 -R "${ref_fasta}" \
    #                 -V ${out_dir}/${Tumor_sample[$i]}/${Tumor_sample[$i]}.vcf \
    #                 -O ${out_dir}/${Tumor_sample[$i]}/${Tumor_sample[$i]}.filter.vcf \
    #                 --filter-expression "QD < 2.0" --filter-name "QD2" \
    #                 --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    #                 --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    #                 --filter-expression "FS > 60.0" --filter-name "FS60" \
    #                 --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    #                 --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    #                 --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

    gatk FilterMutectCalls \
        -V ${out_dir}/${Tumor_sample[$i]}/${Tumor_sample[$i]}.vcf \
        -R ${ref_fasta} \
        -O ${out_dir}/${Tumor_sample[$i]}/${Tumor_sample[$i]}.filter.vcf
done