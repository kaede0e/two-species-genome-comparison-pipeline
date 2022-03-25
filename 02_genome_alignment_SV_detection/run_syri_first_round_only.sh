#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=YOURACCOUNT
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G

genus=$1
echo $genus

#0 Give paths to the command and load modules
export PATH=$PATH:/PATH/TO/anchorwave/minimap2/
export PATH=$PATH:/PATH/TO/bin
export PATH=$PATH:/PATH/TO/bin/syri
export PATH=$PATH:/PATH/TO/YOURDATA

module load scipy-stack
source /PATH/TO/python_env/bin/activate
PATH_TO_SYRI="/PATH/TO/bin/syri/syri/bin/syri"
PATH_TO_PLOTSR="/PATH/TO/bin/syri/syri/bin/plotsr"

python3 $PATH_TO_SYRI \
-c ${genus}_genome/minimap2_ref_qrygenome_chr_alignment.sam \
-r ${genus}_genome/*refgenome_chr.fa \
-q ${genus}_genome/qrygenome_chr_rev.fa \
--dir ${genus}_genome \
-k -F S 2> ${genus}_genome/syri_error_out.txt;
cat ${genus}_genome/syri_error_out.txt | grep WARN  | grep "high fraction of inverted" | cut -f 23 -d " " | sed s/\(//g | sed s/\)\.//g > ${genus}_genome/chr_to_rev.txt ;
echo "Done $genus first SyRI comparison"
