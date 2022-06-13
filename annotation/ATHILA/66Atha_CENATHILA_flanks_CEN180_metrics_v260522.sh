#!/bin/bash

source  activate R-4.0.0
./66Atha_CENATHILA_flanks_CEN180_metrics_v260522.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 1000 5e3 Flanks
conda deactivate
