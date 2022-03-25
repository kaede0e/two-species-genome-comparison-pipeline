#!/bin/bash

## Export path to the command 
PATH_TO_PLOTSR="/PATH/TO/bin/syri/syri/bin/plotsr"

# Visualize the alignment in PDF format (built-in SyRI function)
python3 $PATH_TO_PLOTSR syri.out refgenome_chr.fa qrygenome_chr.fa -H 8 -W 5
