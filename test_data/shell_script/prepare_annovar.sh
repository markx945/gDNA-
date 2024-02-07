#! /bin/bash

## 下载最新版annovar
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz -O ../other_ref/annovar.latest.tar.gz

## 进入目标路径
cd ../other_ref

## 解压缩文件
tar -zxvf annovar.latest.tar.gz

## 进入目标路径，下载注释文件信息
cd annovar

## 下载hg38 refgene注释文件
perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/