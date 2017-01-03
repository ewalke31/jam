#! /bin/bash
# usage: ./bam_variants.bash < path to .sam > < path to ref genoms >
# sample usage: ./bam_variants.bash './tests/ex1' ./tests/ex1.fa

samtools view -b -S -o $1.bam $1.sam
samtools sort $1.bam > $1.sorted
samtools index $1.sorted
samtools mpileup -g -f $2 $1.sorted > $1_variants_bam.bcf
bcftools call -c -v $1_variants_bam.bcf > $1_variants_bam.vcf
