#! /bin/bash

recal_bam=$1

# 首先使用samtools查找高覆盖度区域
samtools depth -a  ${recal_bam} | awk '$3 > 10 {print $1"\t"$2-1"\t"$2}' > high_coverage_10X.bed

# 使用bedtools合并接近的区域
bedtools merge -i high_coverage_10X.bed > merged_10X_coverage.bed