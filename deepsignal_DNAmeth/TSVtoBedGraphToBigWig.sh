#!/bin/bash

# Convert TSV file containing deepsignal-calculated
# context-specific DNA methylation frequencies into bigWig format
# suitable for use with deepTools computeMatrix for fine-scale profiling

# Usage:
# ./TSVtoBedGraphToBigWig.sh Col-0_deepsignalDNAmeth_20kb_MappedOn Col-0.ragtag_scaffolds CpG

acc=$1
refbase=$2
context=$3

source activate BSseq_mapping

#[ -d ./bg ] || mkdir ./bg
#[ -d ./bw ] || mkdir ./bw

# USAGE: bedGraphToBigWig in.bedGraph chrom.sizes out.bw                                                  
# where in.bedGraph is a four-column file in the format:                                                  
#       <chrom> <start> <end> <value>
# and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>                             
# and out.bw is the output indexed big wig file.
# The input bedGraph file must be sorted, use the unix sort command:                                      
(cat ${acc}_${refbase}_${context}.tsv | LC_COLLATE=C sort -k1,1 -k2,2n \
| awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $2, $2+1, $10*100}' - \
> bg/${acc}_${refbase}_${context}.bedGraph; \
bedGraphToBigWig bg/${acc}_${refbase}_${context}.bedGraph /home/ajt200/analysis/pancentromere/assemblies/${refbase}.fa.sizes bw/${acc}_${refbase}_${context}.bw ) \
&> ${acc}_${refbase}_${context}_TSVtoBedGraphToBigWig.log

conda deactivate
