#! /bin/bash


bed_file="/home/cfff_r2636/data/gDNA_231114/high_coverage_10X.bed"
ref_flat="/home/cfff_r2636/data/reference/hg38/refFlat.txt"
tumor_bam="/home/cfff_r2636/data/gDNA_231114/mapping/AF_1_combined/AF_1_combined.recal.bam"
normal_bam="/home/cfff_r2636/data/gDNA_231114/mapping/M8_1_combined/M8_1_combined.recal.bam"
nt="16"
ref_dir="/home/cfff_r2636/data/reference/hg38/genome/hg38.fa"

sample="AF_1_combined"
outdir="/home/cfff_r2636/data/gDNA_231114/CNV"
center="0"

sample2="B34_AF_1"


cnvkit.py target ${bed_file} --annotate ${ref_flat} --split --short-names -o ${outdir}/my_baits.bed
                   
cnvkit.py batch ${tumor_bam} \
                      --normal ${normal_bam} \
                      --targets ${outdir}/my_baits.bed \
                      --fasta ${ref_dir} \
                      --annotate ${ref_flat} -p $nt \
                      --drop-low-coverage \
                      --output-dir ${outdir}/${sample}.reference.cnn

cnvkit.py batch  ${tumor_bam} \
                  -r ${outdir}/${sample}.reference.cnn/reference.cnn \
                  --output-dir ${outdir}/${sample}.cns \
                  -p $nt

# cnvkit.py call ${outdir}/${sample}.cns/${sample}.cns --center-at $center \
#                       -o ${outdir}/${sample}.call.cns

# cnvkit.py scatter ${outdir}/${sample}.cns/${sample2}.cnr -s ${outdir}/${sample2}.call.cns -o ${outdir}/${sample}.scatter.pdf
# cnvkit.py diagram ${outdir}/${sample}.cns/${sample2}.cnr -s ${outdir}/${sample2}.call.cns -o ${outdir}/${sample}.diagram.pdf
# cnvkit.py heatmap ${outdir}/${sample}.cns/${sample2}.cnr ${outdir}/${sample2}.call.cns -o ${outdir}/${sample}.heatmap.pdf