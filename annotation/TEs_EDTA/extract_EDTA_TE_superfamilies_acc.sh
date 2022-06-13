#!/bin/bash

source activate R-4.0.0

./extract_EDTA_TE_superfamilies_acc.R Col-0.ragtag_scaffolds 'Chr1,Chr2,Chr3,Chr4,Chr5' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Col-0.ragtag_scaffolds 'Chr1' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Col-0.ragtag_scaffolds 'Chr2' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Col-0.ragtag_scaffolds 'Chr3' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Col-0.ragtag_scaffolds 'Chr4' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Col-0.ragtag_scaffolds 'Chr5' Gypsy_LTR

./extract_EDTA_TE_superfamilies_acc.R Cvi-0.ragtag_scaffolds 'Chr1,Chr2,Chr3,Chr4,Chr5' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Cvi-0.ragtag_scaffolds 'Chr1' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Cvi-0.ragtag_scaffolds 'Chr2' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Cvi-0.ragtag_scaffolds 'Chr3' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Cvi-0.ragtag_scaffolds 'Chr4' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Cvi-0.ragtag_scaffolds 'Chr5' Gypsy_LTR

./extract_EDTA_TE_superfamilies_acc.R Ler-0_110x.ragtag_scaffolds 'Chr1,Chr2,Chr3,Chr4,Chr5' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Ler-0_110x.ragtag_scaffolds 'Chr1' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Ler-0_110x.ragtag_scaffolds 'Chr2' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Ler-0_110x.ragtag_scaffolds 'Chr3' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Ler-0_110x.ragtag_scaffolds 'Chr4' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Ler-0_110x.ragtag_scaffolds 'Chr5' Gypsy_LTR

./extract_EDTA_TE_superfamilies_acc.R Tanz-1.patch.scaffold.Chr 'Chr1,Chr2,Chr3,Chr4,Chr5' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Tanz-1.patch.scaffold.Chr 'Chr1' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Tanz-1.patch.scaffold.Chr 'Chr2' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Tanz-1.patch.scaffold.Chr 'Chr3' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Tanz-1.patch.scaffold.Chr 'Chr4' Gypsy_LTR
./extract_EDTA_TE_superfamilies_acc.R Tanz-1.patch.scaffold.Chr 'Chr5' Gypsy_LTR

conda deactivate
