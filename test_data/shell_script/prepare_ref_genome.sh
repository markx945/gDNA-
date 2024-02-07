#! /bin/bash

## 下载参考基因组
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O  ../other_ref/hg38.fa.gz

## 下载已知SNP位点

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz -O ../other_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz

## 为参考基因组添加索引
cd ../other_ref

## 解压缩fasta文件
gunzip hg38.fa.gz

## 使用bwa添加索引
bwa index hg38.fa

## 使用samtools添加索引
samtools faidx hg38.fa

echo "Complete genome download and indexing"

