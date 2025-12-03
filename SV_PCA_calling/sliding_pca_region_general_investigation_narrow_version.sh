#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=SPRGI  # sensible name for the job
#SBATCH --mem=50G                 # Default memory per CPU is 3GB. 
#SBATCH --output=SPRGI_bis.txt  #output



echo "
  _____________________________  ________.___ 
 /   _____/\______   \______   \/  _____/|   |
 \_____  \  |     ___/|       _/   \  ___|   |
 /        \ |    |    |    |   \    \_\  \   |
/_______  / |____|    |____|_  /\______  /___|
        \/                   \/        \/     "

echo -e "Perform general steps for interesting regions.\n"



##test parameters

if [ -z "$1" ]
then
      echo "Error: No options provided. Terminating."
      exit
else
      CHROM=$1
      if [ "$CHROM" == "help" ]
      then 
            echo -e "This script is used to perform general investigation of specific regions."
            echo -e "How to use:\n"
            echo -e "sbatch sliding_pca_region_general_investigation_narrow_version.slurm [Chrom] [start_region] [end_region]\n"
            exit
      else
            if [ -z "$2" ]
            then
                  echo "Error: No start provided. Terminating."
                  exit
            else
                  START=$2
                  if [ -z "$3" ]
                  then
                        echo "Error: No end provided. Terminating."
                        exit
                  else
                        END=$3
                  fi
            fi

      fi
fi

WD=data #shortcut for the working directory

DD=datadirectory #shortcut for the data directory

echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."

module load R/4.3.1 

## In this version, the input is the "exact" region we want to investigate

### So first, make vcf of the specified region, and a vcf for the region +- a certain size

# Calculate the result of (END - START)
difference=$((END - START))

# Calculate the size we want to increase the region (here 20%)
size_to_increase=$((difference * 50 / 100))

# Add the size to END
NEWEND=$((END + size_to_increase))

# Subtract the size from START
NEWSTART=$((START - size_to_increase))




## keep usable samples

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -r ${CHROM}:${START}-${END} -S $WD/samples_int_region_bis.txt -O z -o $WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int.vcf.gz $DD/merged_SNP_filtered.recode.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int.vcf.gz


singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -r ${CHROM}:${NEWSTART}-${NEWEND} -S $WD/samples_int_region_bis.txt -O z -o $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz $DD/merged_SNP_filtered.recode.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz


## make also phased vcf

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -r ${CHROM}:${START}-${END} -S $WD/samples_int_region_bis.txt -O z -o $WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int_phased.vcf.gz $DD/subset368_SNP_phased.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int_phased.vcf.gz


singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -r ${CHROM}:${NEWSTART}-${NEWEND} -S $WD/samples_int_region_bis.txt -O z -o $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int_phased.vcf.gz $DD/subset368_SNP_phased.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int_phased.vcf.gz


## extract genotype of the specified regions
sbatch $HOME/get_vcf_stat.slurm $WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int.vcf.gz 258


sbatch $HOME/get_vcf_stat.slurm $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz 258

##make a map file
module load PLINK/1.9b_6.17-x86_64

plink --vcf $WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int.vcf.gz --recode --double-id --aec --out $WD/sprgi_mapfiles/mapfile_${CHROM}_${START}_${END}



## identify inversion breakpoint and cluster individuals on R

## first wait util the stat file exist
file_path=$WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int.vcf.gz_258.stat
timeout_seconds=3000  # Adjust the timeout value as needed

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

echo "File has been generated. Progressing to next step."

echo "Lauching script to identify boundaries and genotype for region ${CHROM}_${START}_${END}"

Rscript $HOME/SPRGI_new_version_getgenotype.R ${CHROM}_${START}_${END}

## wait again until the Rscript is done (can take a long time)

file_path2=$WD/sprgi_gt_AABB/${CHROM}_${START}_${END}_II_ind.txt
timeout_seconds2=10000 # Adjust the timeout value as needed

start_timeb=$(date +%s)

while [ ! -e "$file_path2" ]; do

      echo "Waiting for file to be created...."

    current_time=$(date +%s)
    elapsed_timeb=$((current_time - start_timeb))

    if [ "$elapsed_timeb" -ge "$timeout_seconds2" ]; then
        echo "Timeout reached. File not found. Terminating"
        exit 
    fi

    sleep 800  # Adjust the sleep duration as needed
done

echo "Genotype for AA and BB have been generated. Progressing to next step."


### compute dxy and pi

##make the vcf of invariant + variant site

mkdir $WD/tmpdirbcftools${CHROM}_${NEWSTART}_${NEWEND}

echo "making invariant vcf"

IFD=ifddatadir/ ##directory with a text file of list of all bam files (bamfileslist.txt) used for the software

## use the bam list and the reference genome to make the invariant file
singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools mpileup -f $DD/Simon_Final2021_CHR.fasta -b $IFD/bamfileslist.txt -r ${CHROM}:${NEWSTART}-${NEWEND} | singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools call -m -Oz -f GQ -o $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}.vcf.gz

## subset samples

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -S $IFD/allsample369.txt -r ${CHROM}:${NEWSTART}-${NEWEND} $DD/merged_SNP_filtered.recode.vcf.gz -Oz -o $WD/invariantfilebychrom/merged_SNP_filtered_369_${CHROM}_${NEWSTART}_${NEWEND}.recode.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/invariantfilebychrom/merged_SNP_filtered_369_${CHROM}_${NEWSTART}_${NEWEND}.recode.vcf.gz

##reheader

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools reheader -s $IFD/renaming_file_invariant.txt $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}.vcf.gz -o $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}_renamed.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}_renamed.vcf.gz

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -S $IFD/allsample369.txt $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}_renamed.vcf.gz -Oz -o $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}_renamed_ordered.vcf.gz

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}_renamed_ordered.vcf.gz

## merge invariant and variant vcfs and sort them

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools concat $WD/invariantfilebychrom/invariantvcffile_${CHROM}_${NEWSTART}_${NEWEND}_renamed_ordered.vcf.gz $WD/invariantfilebychrom/merged_SNP_filtered_369_${CHROM}_${NEWSTART}_${NEWEND}.recode.vcf.gz -Oz -o $WD/invariantfilebychrom/merged_SNP_invariant_${CHROM}_${NEWSTART}_${NEWEND}.vcf.gz

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools sort $WD/invariantfilebychrom/merged_SNP_invariant_${CHROM}_${NEWSTART}_${NEWEND}.vcf.gz -O z -o $WD/invariantfilebychrom/merged_SNP_invariant_${CHROM}_${NEWSTART}_${NEWEND}_sorted.vcf.gz --temp-dir $WD/tmpdirbcftools${CHROM}_${NEWSTART}_${NEWEND}

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/invariantfilebychrom/merged_SNP_invariant_${CHROM}_${NEWSTART}_${NEWEND}_sorted.vcf.gz
    


echo "Running Pixy"



## conda activation
echo "Activating Miniconda3 module for $USER"
module load Miniconda3 
eval "$(conda shell.bash hook)"

##run relate to convert vcf to the right format

conda activate pixy

pixy --stats pi dxy \
--vcf $WD/invariantfilebychrom/merged_SNP_invariant_${CHROM}_${NEWSTART}_${NEWEND}_sorted.vcf.gz \
--populations $WD/sprgi_gt_pca_files/gtnopop_${CHROM}_${START}_${END}.txt \
--window_size 10000 \
--output_folder $WD/sprgi_pidxy

mv $WD/sprgi_pidxy/pixy_pi.txt $WD/sprgi_pidxy/pixy_pi_${CHROM}_${NEWSTART}_${NEWEND}.txt
mv $WD/sprgi_pidxy/pixy_dxy.txt $WD/sprgi_pidxy/pixy_dxy_${CHROM}_${NEWSTART}_${NEWEND}.txt




## compute FST between genotypes


echo "Generating FST between genotypes."

singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz --weir-fst-pop $WD/sprgi_gt_AABB/${CHROM}_${START}_${END}_II_ind.txt --weir-fst-pop $WD/sprgi_gt_AABB/${CHROM}_${START}_${END}_NINI_ind.txt --out $WD/sprgi_fstfiles/${CHROM}_${NEWSTART}_${NEWEND}_fst_II_vs_NINI.txt
singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz --weir-fst-pop $WD/sprgi_gt_AABB/${CHROM}_${START}_${END}_II_ind.txt --weir-fst-pop $WD/sprgi_gt_AABB/${CHROM}_${START}_${END}_NII_ind.txt --out $WD/sprgi_fstfiles/${CHROM}_${NEWSTART}_${NEWEND}_fst_II_vs_NII.txt
singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz --weir-fst-pop $WD/sprgi_gt_AABB/${CHROM}_${START}_${END}_NII_ind.txt --weir-fst-pop $WD/sprgi_gt_AABB/${CHROM}_${START}_${END}_NINI_ind.txt --out $WD/sprgi_fstfiles/${CHROM}_${NEWSTART}_${NEWEND}_fst_NII_vs_NINI.txt


## compute FST between populations

echo "Generating FST populations."

singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz --weir-fst-pop $WD/FE_ind.txt --weir-fst-pop $WD/FA_ind.txt --out $WD/sprgi_fstfiles/${CHROM}_${NEWSTART}_${NEWEND}_fst_FE_vs_FA.txt
singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz --weir-fst-pop $WD/WE_ind.txt --weir-fst-pop $WD/WA_ind.txt --out $WD/sprgi_fstfiles/${CHROM}_${NEWSTART}_${NEWEND}_fst_WE_vs_WA.txt



#singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $WD/sprgi_vcffiles/PCA_${CHROM}_${NEWSTART}_${NEWEND}_int.vcf.gz --weir-fst-pop $WD/WE_ind.txt --weir-fst-pop $WD/FE_ind.txt --out $WD/sprgi_fstfiles/${CHROM}_${NEWSTART}_${NEWEND}_fst_FE_vs_WE.txt



### compute PCA on the region
echo "Generating PCA on the region."

module load PLINK/2.00a2.3_x86_64


plink2 --vcf $WD/sprgi_vcffiles/PCA_${CHROM}_${START}_${END}_int.vcf.gz --pca --out $WD/sprgi_pcafiles/PCA_${CHROM}_${START}_${END} --aec --dog 


## plots everything on R



echo "plotting on R"


Rscript $HOME/SPRGI_new_version_plots.R ${CHROM}_${START}_${END}


echo "Script is done"
