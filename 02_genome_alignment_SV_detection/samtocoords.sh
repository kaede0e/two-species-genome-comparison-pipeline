#!/bin/bash

## This should be done on interactive node. 
module load scipy-stack
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/python_env/bin/activate
python3       #the three ">" appears which tells you that you're writing python command
>>> from syri.pyxFiles.synsearchFunctions import samtocoords
>>> table=samtocoords("PATH/TO/SAM")
>>> table.to_csv(r'table_original.txt', index = False, sep = '\t') #makes table.txt file in CSV format in my directory
exit()
