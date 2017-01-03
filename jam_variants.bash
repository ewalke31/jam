#! /bin/bash
# usage: ./jam_variants.bash < path to .sam > < path to ref genoms >
# sample usage: ./jam_variants.bash './tests/ex1' ./tests/ex1.fa

python3 jam.py e v $1.sam $2
python3 jam.py d v $1.vjam $2
samtools view -b -S -o $1.v.bam $1.v.sam
samtools sort $1.v.bam > $1.v.sorted
samtools index $1.v.sorted
samtools mpileup -g -f $2 $1.v.sorted > $1_variants_jam.bcf
bcftools call -c -v $1_variants_jam.bcf > $1_variants_jam.vcf

