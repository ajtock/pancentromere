#!/bin/bash

source activate R-4.0.0

./CEN180_CENAthila_nonCENAthila_Gypsy_DNAmethAllContexts_3metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 180 2000 2000 2000 2kb 10 10 10bp 10bp '0.02,0.96' 'Col-0_deepsignalDNAmeth_20kb_MappedOn_Col-0.ragtag_scaffolds_CpG,Col-0_deepsignalDNAmeth_20kb_MappedOn_Col-0.ragtag_scaffolds_CHG,Col-0_deepsignalDNAmeth_20kb_MappedOn_Col-0.ragtag_scaffolds_CHH' 'pancentromere/deepsignal_DNAmeth/Col-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Col-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Col-0.ragtag_scaffolds' 'Col-0 CG,Col-0 CHG,Col-0 CHH' 'dodgerblue4,dodgerblue1,cyan2' Col-0.ragtag_scaffolds

./CEN180_CENAthila_nonCENAthila_Gypsy_DNAmethAllContexts_3metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 180 2000 2000 2000 2kb 10 10 10bp 10bp '0.02,0.96' 'Cvi-0_deepsignalDNAmeth_10kb_MappedOn_Cvi-0.ragtag_scaffolds_CpG,Cvi-0_deepsignalDNAmeth_10kb_MappedOn_Cvi-0.ragtag_scaffolds_CHG,Cvi-0_deepsignalDNAmeth_10kb_MappedOn_Cvi-0.ragtag_scaffolds_CHH' 'pancentromere/deepsignal_DNAmeth/Cvi-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Cvi-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Cvi-0.ragtag_scaffolds' 'Cvi-0 CG,Cvi-0 CHG,Cvi-0 CHH' 'dodgerblue4,dodgerblue1,cyan2' Cvi-0.ragtag_scaffolds

./CEN180_CENAthila_nonCENAthila_Gypsy_DNAmethAllContexts_3metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 180 2000 2000 2000 2kb 10 10 10bp 10bp '0.02,0.96' 'Ler-0_deepsignalDNAmeth_10kb_MappedOn_Ler-0_110x.ragtag_scaffolds_CpG,Ler-0_deepsignalDNAmeth_10kb_MappedOn_Ler-0_110x.ragtag_scaffolds_CHG,Ler-0_deepsignalDNAmeth_10kb_MappedOn_Ler-0_110x.ragtag_scaffolds_CHH' 'pancentromere/deepsignal_DNAmeth/Ler-0_110x.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Ler-0_110x.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Ler-0_110x.ragtag_scaffolds' 'Ler-0 CG,Ler-0 CHG,Ler-0 CHH' 'dodgerblue4,dodgerblue1,cyan2' Ler-0_110x.ragtag_scaffolds

./CEN180_CENAthila_nonCENAthila_Gypsy_DNAmethAllContexts_3metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 180 2000 2000 2000 2kb 10 10 10bp 10bp '0.02,0.96' 'Tanz-1_deepsignalDNAmeth_unkb_MappedOn_Tanz-1.patch.scaffold.Chr_CpG,Tanz-1_deepsignalDNAmeth_unkb_MappedOn_Tanz-1.patch.scaffold.Chr_CHG,Tanz-1_deepsignalDNAmeth_unkb_MappedOn_Tanz-1.patch.scaffold.Chr_CHH' 'pancentromere/deepsignal_DNAmeth/Tanz-1.patch.scaffold.Chr,pancentromere/deepsignal_DNAmeth/Tanz-1.patch.scaffold.Chr,pancentromere/deepsignal_DNAmeth/Tanz-1.patch.scaffold.Chr' 'Tanz-1 CG,Tanz-1 CHG,Tanz-1 CHH' 'dodgerblue4,dodgerblue1,cyan2' Tanz-1.patch.scaffold.Chr

conda deactivate
