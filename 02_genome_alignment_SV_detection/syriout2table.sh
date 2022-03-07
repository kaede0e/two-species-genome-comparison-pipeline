#!/bin/bash

#Format SyRI output into tidy table
cat invOut.txt | perl /PATH/TO/syriout2table.pl > invOut_table.txt
cat synOut.txt | perl /PATH/TO/syriout2table.pl > synOut_table.txt
cat TLOut.txt | perl /PATH/TO/syriout2table.pl > TLOut_table.txt
cat invTLOut.txt | perl /PATH/TO/syriout2table.pl > invTLOut_table.txt
cat invDupOut.txt | perl /PATH/TO/syriout2table.pl > invDupOut_table.txt
cat dupOut.txt | perl /PATH/TO/syriout2table.pl > dupOut_table.txt
