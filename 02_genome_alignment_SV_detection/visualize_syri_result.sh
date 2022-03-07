#!/bin/bash

## Export path to the command 
PATH_TO_PLOTSR="/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/syri/syri/bin/plotsr"

# Visualize the alignment in PDF format (built-in SyRI function)
python3 $PATH_TO_PLOTSR syri.out refgenome_chr.fa  -H 8 -W 5
