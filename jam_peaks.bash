#! /bin/bash
# usage: ./jam_variants.bash < path to .sam > < path to ref genoms >
# sample usage: ./jam_variants.bash './tests/ex1' ./tests/ex1.fa

python3 jam.py e p $1.sam $2
python3 jam.py d p $1.pjam $2
samtools view -bS -o $1.p.bam $1.p.sam
MACS14/bin/macs14 -t $1.p.bam --nomodel -g dm -w --space=20 -S -n $1.p
