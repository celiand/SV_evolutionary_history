#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=fstclus  # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=fstclus.txt  #output


WD=data/ #shortcut for the working directory

echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."

DD=directorydata/ ### data with a vcf containing all SNPs for all individuals

TP=data/tempfilefst

RD=data/fst_files_cluster

module load R/4.3.1 

## get all unique code 

CODES=($(tail -n +2 $HOME/full_info_windows_PCA.txt | awk '{print $9}' | sort -u))


## perform fst for each code value (each region)

for value in ${CODES[@]}; do
        ## launch Rscript to extract individuals belonging to each group and region coordinates

        Rscript data/fst_cluster_process.R $value

        value=${value//\"/}

        ##check if results are done (coordinates files present)
        var=true

    while $var ; do
        if ls ${TP}/*coordinates.txt >/dev/null 2>&1; then
            var=false
            echo "Files found in the directory. Continuing process"
            break

        else
            echo "Files not found, waiting for rscript to finish"
            sleep 200  # Sleep for 200 seconds before rechecking
        fi
    done

    file_count=$(ls -l "$TP" | grep "^-" | wc -l)

    ## substract 1 to not take into account the coordinate file 
    file_count=$((file_count - 1))

    ### substract vcf 
    echo "substracting vcf"
    singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -R $TP/${value}_coordinates.txt -O z -o $TP/temp_region.vcf.gz $DD/merged_SNP_filtered.recode.vcf.gz

    singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p -f vcf $TP/temp_region.vcf.gz

    echo "doing fst for each combination"
    for ((i = 1; i <= $file_count; i++)); do
        for ((j = i + 1; j <= $file_count; j++)); do
            ## perform fst for each combinations
            #echo "command -option $i -option $j"
            singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $TP/temp_region.vcf.gz --weir-fst-pop $TP/${value}_${i}.txt --weir-fst-pop $TP/${value}_${j}.txt --out $RD/FST_${value}_${i}_${j}
        done
    done
    cd $TP
    rm *.txt
    rm *.vcf*
    cd $WD


done



