#!/bin/bash

#SBATCH --cpus-per-task=2            # 1 core(CPU)
#SBATCH --job-name=formatslim  # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=formatslim.txt  #output


WD=Data/slimfolder


echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."


echo -e "POS\tMutation_type\tFrequency\tPopulation\tRun\tType_simu">$WD/all_vcf_info.txt

#### list all type of simulations
Arraysimu=("1.1_500_recurrent_sv_selec" "1.1_1000_recurrent_sv_selec" "1.1_1500_recurrent_sv_selec" "1.1_500_recurrent_sv_selec_het" "1.1_1000_recurrent_sv_selec_het" "1.1_1500_recurrent_sv_selec_het" "0.5_1500_recurrent_sv_selec_het" "1.1_ancient_old_sv_selec_het" "1.1_ancient_old_sv_selec")

# Define the index you want to extract
CurrentType=$(( SLURM_ARRAY_TASK_ID - 1 )) ##do -1 because array start at 0

# Extract the Xth value
typesimu=${Arraysimu[$CurrentType]}

module load VCFtools

for file in $WD/vcf_fromslim/*${typesimu}.vcf; do
    if [ -f "$file" ]; then
        echo "Processing $file"
        pop_number=$(echo "$file" | sed -n 's/.*pop\([0-9]\+\).*/\1/p')

        run_number=$(echo "$file" | sed -n 's/.*run\([0-9]\+\).*/\1/p')

        echo "Extracting vcf info"
        bcftools query -f "%POS %INFO/MT\n" $file -o temp_vcf${typesimu}.stat

        echo "Extracting variant frequency"
        vcftools --vcf $file --freq --out temp_vcf${typesimu}.vaf

        nbofSNPs=$(wc -l temp_vcf${typesimu}.stat | awk '{print $1}')
        ##echo $nbofSNPs
        for i in $(seq 1 $nbofSNPs); do

            posSNP=$(awk -v i="$i" 'NR==i {print $1}' temp_vcf${typesimu}.stat)
            MTSNP=$(awk -v i="$i" 'NR==i {print $2}' temp_vcf${typesimu}.stat)

            FreqSNP=$(awk -v i="$posSNP" '$2==i {print $6}' temp_vcf${typesimu}.vaf.frq | cut -d':' -f2)

            echo -e "$posSNP\t$MTSNP\t$FreqSNP\t$pop_number\t$run_number\t$typesimu">>$WD/all_vcf_info.txt
        done

    fi
done