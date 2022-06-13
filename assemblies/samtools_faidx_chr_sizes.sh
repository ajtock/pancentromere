#!/bin/bash

for i in *.fa; do
#  /home/ajt200/miniconda3/envs/ChIPseq_mapping/bin/samtools faidx ${i}
#cut -f1,2 ${i}.fai > ${i}.chrom.sizes
ln -s ${i}.chrom.sizes ${i}.sizes 
done
