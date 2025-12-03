#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=manta   # sensible name for the job
#SBATCH --mem=35G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=manta_job_output.%A_%a.txt  #output
#SBATCH --partition=smallmem,hugemem   ## work on whatever node available

WD="my_working_directory" #shortcut for the working directory
TMPWD="my_temp_directory" ### set up a temporary directory

REGION=${1} #### region from which the sample is from
$SAMPLENAME=${2} #### name of the sample


### make a chromosome file name (different chromosome naming convention depending on the data)

	if [ "$REGION" == "Wild_European" ] ##rename bam file if needed
	then
		cp /genomedir/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz $TMPWD/NCBI_CHR.fasta
		cp /genomedir/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz.fai $TMPWD/NCBI_CHR.fasta.fai
		cp /genomedir/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz.gzi $TMPWD/NCBI_CHR.fasta.gzi
	  		echo "#CHROM	START	END">$TMPWD/chromregion.bed
			  	for ((i=42; i<71; i++))
				do
					CHROM="NC_0594"$i".1"
					echo "$CHROM 	1	300000000">>$TMPWD/chromregion.bed 
				done
	  		
	  		singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif bgzip -f $TMPWD/chromregion.bed
	  		singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p bed $TMPWD/chromregion.bed.gz
	elif [ "$REGION" == "Wild_American" ]
	then
		#ref genome and index
		cp /genomedir/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz $TMPWD/NCBI_CHR.fasta
		cp /genomedir/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz.fai $TMPWD/NCBI_CHR.fasta.fai
		cp /genomedir/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz.gzi $TMPWD/NCBI_CHR.fasta.gzi

		echo "#CHROM	START	END">$TMPWD/chromregion.bed
			  	for ((i=42; i<71; i++))
				do
					CHROM="NC_0594"$i".1"
					echo "$CHROM 	1	300000000">>$TMPWD/chromregion.bed 
				done
	  		
	  		singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif bgzip -f $TMPWD/chromregion.bed
	  		singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p bed $TMPWD/chromregion.bed.gz
	else
		#ref genome and index
		cp $WD/data/Simon_Final2021_CHR.fasta $TMPWD/Simon_Final2021_CHR.fasta
		cp $WD/data/Simon_Final2021_CHR.fasta.fai $TMPWD/Simon_Final2021_CHR.fasta.fai
	fi

	



### 1 by 1 analysis ###

## run the config manta command (not the same reference genome depending on chromosome naming)

	
	if [ "$REGION" == "Wild_European" ] ##rename bam file if needed
	then
		singularity exec $TMPWD/manta_1_6_2_h9ee0642_2.sif configManta.py --bam $TMPWD/$SAMPLENAME.bam --generateEvidenceBam --referenceFasta $TMPWD/NCBI_CHR.fasta --callRegions $TMPWD/chromregion.bed.gz --runDir $TMPWD 
	elif [ "$REGION" == "Wild_American" ]
	then
		singularity exec $TMPWD/manta_1_6_2_h9ee0642_2.sif configManta.py --bam $TMPWD/$SAMPLENAME.bam --generateEvidenceBam --referenceFasta $TMPWD/NCBI_CHR.fasta --callRegions $TMPWD/chromregion.bed.gz --runDir $TMPWD 
	else
		singularity exec $TMPWD/manta_1_6_2_h9ee0642_2.sif configManta.py --bam $TMPWD/$SAMPLENAME.bam --generateEvidenceBam --referenceFasta $TMPWD/Simon_Final2021_CHR.fasta --runDir $TMPWD 
	fi
	
    ### run sv detection
	singularity exec $TMPWD/manta_1_6_2_h9ee0642_2.sif $TMPWD/runWorkflow.py

## save results
cp $TMPWD/results/variants/diploidSV.vcf.gz /storage/vcfcall/manta_${ANALYSIS}_${REGION}_${SAMPLENAME}.vcf.gz
cp $TMPWD/results/variants/diploidSV.vcf.gz.tbi /storage/vcfcall/manta_${ANALYSIS}_${REGION}_${SAMPLENAME}.vcf.gz.tbi


	
