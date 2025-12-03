#!/bin/bash

#SBATCH --cpus-per-task=2            # 1 core(CPU)
#SBATCH --job-name=ibshaplotype  # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. 
#SBATCH --output=ibshaplotype.txt  #output


WD=Data/slimfolder


echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."

## transfer the files from "sampledkept.txt" (see analyze_svsimu.R) to this directory
resdir="$WD/vcfkeptanalysis"

for pop in 2 3 4 5; do
  awk -v pop=$pop -v resdir="$resdir" '
    NR>1 {
      run=$2; type=$1;
      fname="pop" pop "_vcf_run" run "_" type ".vcf";
      print "cp vcf_fromslim/" fname " " resdir "/";
    }' sampledkept.txt
done | bash

fi

#### commands to make the individuals sample name:
##  bcftools query -l pop2_vcf_run1286_1.1_1500_recurrent_sv_selec.vcf.gz > raw_sample_simu.txt
## sed 's/$/_pop2/' raw_sample_simu.txt > pop2_sample_simu.txt

combinations=$(for file in "$resdir"/*.vcf; do
    filename="${file%.vcf}"

    # Extract run number (between 'run' and next '_')
    run_number=$(echo "$filename" | sed -E 's/.*run([0-9]+)_.*/\1/')

    # Extract type_simu (everything after runX_)
    type_simu=$(echo "$filename" | sed -E 's/.*run[0-9]+_(.*)/\1/')

    echo "${run_number}|${type_simu}"
done | sort -u)


# Now loop once per unique combination
while IFS="|" read -r run_number type_simu; do
    echo "Processing run=$run_number type=$type_simu"

   

    ### then rename vcf to not have the same individuals name

    bcftools reheader -s pop2_sample_simu.txt $resdir/pop2_vcf_run${run_number}_${type_simu}.vcf -o $resdir/pop2_vcf_run${run_number}_${type_simu}_headerrename.vcf

    bgzip $resdir/pop2_vcf_run${run_number}_${type_simu}_headerrename.vcf

    tabix -p vcf $resdir/pop2_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz

    bcftools reheader -s pop3_sample_simu.txt $resdir/pop3_vcf_run${run_number}_${type_simu}.vcf -o $resdir/pop3_vcf_run${run_number}_${type_simu}_headerrename.vcf

    bgzip $resdir/pop3_vcf_run${run_number}_${type_simu}_headerrename.vcf

    tabix -p vcf $resdir/pop3_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz

    bcftools reheader -s pop4_sample_simu.txt $resdir/pop4_vcf_run${run_number}_${type_simu}.vcf -o $resdir/pop4_vcf_run${run_number}_${type_simu}_headerrename.vcf

    bgzip $resdir/pop4_vcf_run${run_number}_${type_simu}_headerrename.vcf

    tabix -p vcf $resdir/pop4_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz

    bcftools reheader -s pop5_sample_simu.txt $resdir/pop5_vcf_run${run_number}_${type_simu}.vcf -o $resdir/pop5_vcf_run${run_number}_${type_simu}_headerrename.vcf

    bgzip $resdir/pop5_vcf_run${run_number}_${type_simu}_headerrename.vcf

    tabix -p vcf $resdir/pop5_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz


    ### then merge vcf into one

    bcftools merge -o $resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz -O z $resdir/pop5_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz $resdir/pop4_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz $resdir/pop3_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz $resdir/pop2_vcf_run${run_number}_${type_simu}_headerrename.vcf.gz
    
    tabix -p vcf $resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz 


    ### then compute the distance haplotype metric
    echo "haplotype stat"

    module load PLINK/2.00a2.3_x86_64

    # We use plink2 to read the VCF and generate PLINK binary files (.bed/.bim/.fam).
    plink2 \
    --vcf $resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz  \
    --double-id \
    --set-all-var-ids @:#:\$r:\$a \
    --max-alleles 2 \
    --allow-extra-chr \
    --make-bed \
    --out region

    module load PLINK/1.9b_6.17-x86_64

    # Here we use the “classic” plink (v1.9).
    # --distance square 1-ibs produces a square matrix of pairwise distances.
    # “1-IBS” means the proportion of sites that are *not identical-by-state*.
    # The output is:
    #   region_ibs.mdist     : distance matrix
    #   region_ibs.mdist.id  : sample IDs (order of rows/columns)
    plink \
    --bfile region \
    --distance square 1-ibs \
    --out region_ibs

   
    bcftools query -r 1:40001-40001 -f '[%SAMPLE\t%GT\n]' $resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz   | awk '$2=="1/1" || $2=="1|1" {print $1}' > group1_ids.txt

    bcftools query -r 1:40001-40001 -f '[%SAMPLE\t%GT\n]' $resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz   | awk '$2=="0/0" || $2=="0|0" {print $1}' > group2_ids.txt


## get genotype data

    echo "get genotype file"
        sbatch $HOME/get_vcf_stat.slurm $resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz 258

        file_path=$resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz_258.stat
        timeout_seconds=5000  # Adjust the timeout value as needed

        start_time=$(date +%s)

        while [ ! -e "$file_path" ]; do

            echo "Waiting for file to be created...."

            current_time=$(date +%s)
            elapsed_time=$((current_time - start_time))

            if [ "$elapsed_time" -ge "$timeout_seconds" ]; then
                echo "Timeout reached. File not found. Terminating"
                exit 
            fi

            sleep 500  # Adjust the sleep duration as needed
        done

        mv $resdir/allpop_vcf_run${run_number}_${type_simu}.vcf.gz_258.stat $resdir/allpop_vcf_258.stat

        ### Rscript

        echo "rscript to make files"
        Rscript $WD/run_ibsstat.R ${run_number}_${type_simu}


done <<< "$combinations"