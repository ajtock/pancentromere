#!/bin/bash

csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Col-0_deepsignalDNAmeth_20kb_MappedOn Col-0.ragtag_scaffolds CpG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Col-0_deepsignalDNAmeth_20kb_MappedOn Col-0.ragtag_scaffolds CHG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Col-0_deepsignalDNAmeth_20kb_MappedOn Col-0.ragtag_scaffolds CHH" & sleep 10;

csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Cvi-0_deepsignalDNAmeth_10kb_MappedOn Cvi-0.ragtag_scaffolds CpG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Cvi-0_deepsignalDNAmeth_10kb_MappedOn Cvi-0.ragtag_scaffolds CHG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Cvi-0_deepsignalDNAmeth_10kb_MappedOn Cvi-0.ragtag_scaffolds CHH" & sleep 10;

csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Ler-0_deepsignalDNAmeth_10kb_MappedOn Ler-0_110x.ragtag_scaffolds CpG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Ler-0_deepsignalDNAmeth_10kb_MappedOn Ler-0_110x.ragtag_scaffolds CHG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Ler-0_deepsignalDNAmeth_10kb_MappedOn Ler-0_110x.ragtag_scaffolds CHH" & sleep 10;

csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Tanz-1_deepsignalDNAmeth_unkb_MappedOn Tanz-1.patch.scaffold.Chr CpG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Tanz-1_deepsignalDNAmeth_unkb_MappedOn Tanz-1.patch.scaffold.Chr CHG" & sleep 10;
csmit -m 10G -c 1 "bash ./TSVtoBedGraphToBigWig.sh Tanz-1_deepsignalDNAmeth_unkb_MappedOn Tanz-1.patch.scaffold.Chr CHH" & sleep 10;
