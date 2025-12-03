#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=filter_and_merge   # sensible name for the job
#SBATCH --mem=10G                 # Default memory per CPU is 3GB. 
#SBATCH --output=fandmerge.txt  #output

DD="my_data_directory" #shortcut for the working directory
WD="my_working_directory" 
TMPWD="my_temp_directory" ### set up a temporary directory

##uncompress fasta in TMPWD

FILELIST=$(ls /storage/vcfcall/manta_Single_*.vcf.gz)

for vcf in $(ls $FILELIST)
	do
        TMPFILENAME=$(echo ${vcf##*/}) ##cut the first part of the name (directory)
        TMPFILENAME2=$(echo ${TMPFILENAME%%.*}) ##cut the last part of the name
		if echo "$TMPFILENAME" | grep -q "Wild"
		then
			## rename wild samples chromosomes 
			singularity exec containers/bcftools_1_9_ha228f0b_4.sif bcftools query -f '%CHROM\n' /storage/vcfcall/manta_Single_Wild_European_Aroy_15_0180_RG_Wild_European_renamed.vcf.gz | uniq > wild_chromosome.txt
			singularity exec containers/bcftools_1_9_ha228f0b_4.sif bcftools query -f '%CHROM\n' /storage/vcfcall/manta_Single_Europe_2014G_NO_Males_3000_D02_RG_Europe_renamed.vcf.gz | uniq > farmed_chromosome.txt 
			paste wild_chromosome.txt farmed_chromosome.txt > renamed_chromosome.txt
			singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools annotate --rename-chrs $WD/renamed_chromosome.txt $vcf -O v -o TMPWD/$TMPFILENAME2.vcf 
		else
			singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif bgzip -c -d $vcf > TMPWD/$TMPFILENAME2.vcf
		fi
	done

ls $TMPWD/*.vcf > $TMPWD/listsample.txt

##merge all the samples
echo "Merging all samples."
singularity exec $WD/containers/survivor_1_07_he513fc3_0.sif SURVIVOR merge $TMPWD/listsample.txt 100 1 0 0 0 30 $TMPWD/ALL_SAMPLES.vcf


##sort the vcf
echo "Sorting the vcf."
singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools sort $TMPWD/ALL_SAMPLES.vcf -T $TMPWD/tmpbcftools/ -O z -o $TMPWD/ALL_SAMPLES_sorted.vcf.gz

## remove sites with strange name that cause issues downstream
echo "Removing sites with non usual ALT allele name."
singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -H $TMPWD/ALL_SAMPLES_sorted.vcf.gz | awk '$5 !~/[]\[]/ {print $0}' > $TMPWD/ALL_SAMPLES_sorted_trimmed_data.vcf

singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -h $TMPWD/ALL_SAMPLES_sorted.vcf.gz > $TMPWD/ALL_SAMPLES_header.vcf

cat $TMPWD/ALL_SAMPLES_header.vcf $TMPWD/ALL_SAMPLES_sorted_trimmed_data.vcf > $TMPWD/ALL_SAMPLES_sorted_trimmed_full.vcf


### ilter vcf
singularity exec $WD/containers/vcftools_0_1_16_p15321hd03093a_7.sif vcftools --gzvcf $TMPWD/ALL_SAMPLES_sorted_trimmed_full.vcf --minQ 100 --recode --recode-INFO-all --out $TMPWD/ALL_SAMPLES_sorted_trimmed_full_filtered.vcf


## save results
cp $TMPWD/ALL_SAMPLES_sorted_trimmed_full_filtered.vcf /storage/mergevcf

