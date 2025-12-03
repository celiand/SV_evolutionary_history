#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=repeatsv  # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=repeatsv.txt  #output


WD=Data #shortcut for the working directory

echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."


module load R/4.3.1 


for ((i=1; i<1001; i++))
    do

    Rscript $Home/makerandomsv.R ### R script to make random SVs (see SV_TE_overlap.R, first section)

    file_path=Data/random_SV_boundaries_table.bed

    sleep 2

    while [ ! -e "$file_path" ]; do
     
     sleep 5

    done

### overlap the SV with the repeats
### with real data (run once)	
    ##singularity exec containers/bedtools_2_30_hc088bd4_0.sif bedtools intersect -a $WD/SV_boundaries_table.bed -b $HOME/maskeddatafull.bed -wa -wb > SV_boundaries_overlap_repeat.txt
	
##with the random SVs
    singularity exec containers/bedtools_2_30_hc088bd4_0.sif bedtools intersect -a $WD/random_SV_boundaries_table.bed -b $HOME/maskeddatafull.bed -wa -wb > random_SV_boundaries_table_overlap_repeat.txt

    rm $WD/random_SV_boundaries_table.bed

    Rscript $HOME/computerepeatsvstat.R $i ### call a script to compute the porportion of Sv overlapped (see SV_TE_overlap.R)
	
##just wait that the script is done before moving to the next occurence
    file_path2=Data/tempendoccu.txt

    sleep 2

    while [ ! -e "$file_path2" ]; do
     
     sleep 5

    done

    rm Data/tempendoccu.txt

    done





