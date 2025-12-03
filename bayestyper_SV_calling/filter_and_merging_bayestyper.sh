#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=postprocessSV   # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=postprocessSV.txt  #output


WD="my_working_directory"  #shortcut for the working directory

BD=$WD/bayertyper 

# index and compress all vcfs

echo "Compressing vcf"
for i in {1..74}; do
    singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif bgzip -f $BD/genotyped_file_${i}.vcf
    singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $BD/genotyped_file_${i}.vcf.gz
done


echo "Merging vcf"

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools merge --info-rules ACP:max -O z -o bayestyper_ALL_SV_non_filtered.vcf.gz $BD/genotyped_file_1.vcf.gz $BD/genotyped_file_2.vcf.gz $BD/genotyped_file_3.vcf.gz \
 $BD/genotyped_file_4.vcf.gz $BD/genotyped_file_5.vcf.gz $BD/genotyped_file_6.vcf.gz $BD/genotyped_file_7.vcf.gz $BD/genotyped_file_8.vcf.gz $BD/genotyped_file_9.vcf.gz \
 $BD/genotyped_file_10.vcf.gz $BD/genotyped_file_11.vcf.gz $BD/genotyped_file_12.vcf.gz $BD/genotyped_file_13.vcf.gz $BD/genotyped_file_14.vcf.gz $BD/genotyped_file_15.vcf.gz \
 $BD/genotyped_file_16.vcf.gz $BD/genotyped_file_17.vcf.gz $BD/genotyped_file_18.vcf.gz $BD/genotyped_file_19.vcf.gz $BD/genotyped_file_20.vcf.gz $BD/genotyped_file_21.vcf.gz \
 $BD/genotyped_file_22.vcf.gz $BD/genotyped_file_23.vcf.gz $BD/genotyped_file_24.vcf.gz $BD/genotyped_file_25.vcf.gz $BD/genotyped_file_26.vcf.gz $BD/genotyped_file_27.vcf.gz $BD/genotyped_file_28.vcf.gz \
 $BD/genotyped_file_29.vcf.gz $BD/genotyped_file_30.vcf.gz $BD/genotyped_file_31.vcf.gz $BD/genotyped_file_32.vcf.gz $BD/genotyped_file_33.vcf.gz $BD/genotyped_file_34.vcf.gz \
 $BD/genotyped_file_35.vcf.gz $BD/genotyped_file_36.vcf.gz $BD/genotyped_file_37.vcf.gz $BD/genotyped_file_38.vcf.gz $BD/genotyped_file_39.vcf.gz $BD/genotyped_file_40.vcf.gz \
 $BD/genotyped_file_41.vcf.gz $BD/genotyped_file_42.vcf.gz $BD/genotyped_file_43.vcf.gz $BD/genotyped_file_44.vcf.gz $BD/genotyped_file_45.vcf.gz $BD/genotyped_file_46.vcf.gz $BD/genotyped_file_47.vcf.gz \
 $BD/genotyped_file_48.vcf.gz $BD/genotyped_file_49.vcf.gz $BD/genotyped_file_50.vcf.gz $BD/genotyped_file_51.vcf.gz $BD/genotyped_file_52.vcf.gz $BD/genotyped_file_53.vcf.gz \
 $BD/genotyped_file_54.vcf.gz $BD/genotyped_file_55.vcf.gz $BD/genotyped_file_56.vcf.gz $BD/genotyped_file_57.vcf.gz $BD/genotyped_file_58.vcf.gz $BD/genotyped_file_59.vcf.gz \
 $BD/genotyped_file_60.vcf.gz $BD/genotyped_file_61.vcf.gz $BD/genotyped_file_62.vcf.gz $BD/genotyped_file_63.vcf.gz $BD/genotyped_file_64.vcf.gz $BD/genotyped_file_65.vcf.gz $BD/genotyped_file_66.vcf.gz \
 $BD/genotyped_file_67.vcf.gz $BD/genotyped_file_68.vcf.gz $BD/genotyped_file_69.vcf.gz $BD/genotyped_file_70.vcf.gz $BD/genotyped_file_71.vcf.gz $BD/genotyped_file_72.vcf.gz \
 $BD/genotyped_file_73.vcf.gz $BD/genotyped_file_74.vcf.gz 



 ## filter the obtained vcf
echo "Filtering vcf"
singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf bayestyper_ALL_SV_non_filtered.vcf.gz --minQ 10 --max-missing 0.25 --recode --out bayestyper_ALL_SV_filtered

## compress and index
singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif bgzip -f bayestyper_ALL_SV_filtered.recode.vcf
singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf bayestyper_ALL_SV_filtered.recode.vcf.gz
