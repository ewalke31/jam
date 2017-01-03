#! /bin/bash
# usage: ./jam_variants.bash < path to .sam > < path to ref genoms >
# sample usage: ./jam_variants.bash './tests/ex1' ./tests/ex1.fa

samtools view -bS -o $1.bam $1.sam
MACS14/bin/macs14 -t $1.bam --nomodel -g dm -w --space=20 -S -n $1.orig
