#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=bayestyper   # sensible name for the job
#SBATCH --mem=80G                 # Default memory per CPU is 3GB. 
#SBATCH --output=bayestyper.%A_%a.txt  #output
#SBATCH --exclude=cn-5


##https://github.com/bioinformatics-centre/BayesTyper

WD="my_working_directory"  #shortcut for the working directory
DD=/storage/vcfcall/ ##data directory

## Set up which part of the code is runnig or not
STEP0=false
STEP1=false
STEP2=false


# Check if the variable is true
if [ "$STEP0" = true ]; then

    echo "Creating sample file"

## step 0: make a list of all samples 

FILES_FE=($(ls $DD/Farmed_European/*.vcf.gz))
FILES_FA=($(ls $DD/Farmed_American/*.vcf.gz))
FILES_WA=($(ls $DD/Wild_American/*.vcf.gz))
FILES_WE=($(ls $DD/Wild_European/*.vcf.gz))

for file in "${FILES_FE[@]}"
do
  # Extract the filename without the directory path
  filename=$(basename "$file")
  # Remove everything after the last "." character)
  filename=${filename%%.*}  ##cut everything after the dots
  filename=${filename#*_}  ##cut the first thing before the underscore
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename%_*} ##cut the first thing after the underscore
  filename=${filename%_*}
  # Add the modified filename to the new array
  echo "$filename">>$WD/bayertyper/sample_name.txt
done


for file in "${FILES_FA[@]}"
do
  # Extract the filename without the directory path
  filename=$(basename "$file")
  # Remove everything after the last "." character
  filename=${filename%%.*}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename%_*}
  filename=${filename%_*}
  # Add the modified filename to the new array
  echo "$filename">>$WD/bayertyper/sample_name.txt
done


for file in "${FILES_WA[@]}"
do
  # Extract the filename without the directory path
  filename=$(basename "$file")
  # Remove everything after the last "." character
  filename=${filename%%.*}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename%%_*}
  # Add the modified filename to the new array
  echo "$filename">>$WD/bayertyper/sample_name.txt
done

for file in "${FILES_WE[@]}"
do
  # Extract the filename without the directory path
  filename=$(basename "$file")
  # Remove everything after the last "." character
  filename=${filename%%.*}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename#*_}
  filename=${filename%_*}
  filename=${filename%_*}
  filename=${filename%_*}
  # Add the modified filename to the new array
  echo "$filename">>$WD/bayertyper/sample_name.txt
done

echo "sample file is created"

##make chr ploidy file
vcf_file="$DD/Farmed_European/manta_Single_Europe_2017NOS1CBR_7121_H02_RG_Europe_renamed.vcf.gz"
output_file="chr_ploidy.txt" 

# Use bcftools to extract the chromosome/contig names from the header of the VCF file
singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools view -h "$vcf_file" | grep -oP '(?<=##contig=<ID=)[^,]+' > "$output_file"

awk -v OFS='\t' '{print $0, 2, 2}' "$output_file" > tmpfile && mv tmpfile "$output_file"

fi ##end of the step 0


# Check if the variable is true
if [ "$STEP1" = true ]; then

    echo "Generating variants candidates for all samples"

##step 1 : generate a set of variants candidates for all samples

## create tmp dir
TMPWD="my_temp_directory" ### set up a temporary directory


FILES=($( cat $WD/bayertyper/sample_name.txt))

for file in "${FILES[@]}"
do
  echo "file is $file"
  # Remove everything after the last "." character
  filename=$(echo "$file" | sed -e 's/[^[:print:]]//g')
  # Add the modified filename to the new array
  FILENAMES+=("$filename")
done



cd $DD
##loop on all the files
for file in "${FILENAMES[@]}"
do
    fileprocess=$(find -name manta*${file}*vcf.gz)
    ##if it's a wild sample, rename chromosome
    if [[ $fileprocess == *"Wild"* ]]; then
        echo "Wild sample identified."
        singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools annotate --rename-chrs $WD/renamed_chromosome.txt $fileprocess -O z -o $TMPWD/bt_${file}.vcf.gz 
    else
      cp $fileprocess $TMPWD/bt_${file}.vcf.gz
    fi
    singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $TMPWD/bt_${file}.vcf.gz

    ## Convert allele from output
    singularity exec $WD/containers/bayestyper_1_5_h176a8bc_0.sif bayesTyperTools convertAllele -v $TMPWD/bt_${file}.vcf.gz -g $WD/data/Simon_Final2021_CHR.fasta -o $TMPWD/Converted_${file}

    ## normalize output
    singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools norm -f $WD/data/Simon_Final2021_CHR.fasta -o $TMPWD/Converted_normalized_${file}.vcf.gz -O z $TMPWD/Converted_${file}.vcf

    ##get sample name as in the vcf
    samp=$(singularity exec $WD/containers/bcftools_1_9_ha228f0b_4.sif bcftools query -l $TMPWD/bt_${file}.vcf.gz)

    fname="$TMPWD/Converted_normalized_${file}.vcf.gz"
    # Initialize an empty string
    result=""
    # Concatenate the name and file for the current index, separated by a comma
    result+="${samp}:${fname},"

done

# Remove the trailing comma
result="${result%,}"

## combine output
echo "Combining output."
singularity exec $WD/containers/bayestyper_1_5_h176a8bc_0.sif bayesTyperTools combine -o $WD/bayertyper/combined_files -v $result

## compress output and index

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif bgzip -f $WD/bayertyper/combined_files.vcf

singularity exec $WD/containers/tabix_1_11_hdfd78af_0.sif tabix -p vcf $WD/bayertyper/combined_files.vcf.gz


fi ##end of the step 1

##start of step 2 (genotyping)
if [ "$STEP2" = true ]; then

##choose the sample to run on
ALLFILE=($(cat $WD/bayertyper/sample_name.txt))
##369 samples --> 37 batchs of 10 // 74 batchs of 5


start=$((($SLURM_ARRAY_TASK_ID - 1) * 5))
end=$(($SLURM_ARRAY_TASK_ID * 5 - 1))
#end=5 ##structure of extraction below is start:length

# Create a new array to store the selected values
selected_values=("${ALLFILE[@]:$start:$end}")

## create tmp dir
TMPWD="my_temp_directory" ### set up a temporary directory

  ##copy all the files

  cp $WD/containers/bayestyper_1_5_h176a8bc_0.sif ./
  cp $WD/containers/kmc_3_2_1_h95f258a_1.sif ./
  cp $WD/containers/tabix_1_11_hdfd78af_0.sif ./

for file in "${selected_values[@]}"
do

## download the fasta files or copy from storage

  if [[ $file == *"NO"* ]]; then
      ##copy files for storage
      cp mystorage/fastq/*${file}*.fastq.gz ./
  elif [[ $file == *"SRR"* ]]; then
      ##dl file
    #Download fastQ files from NCBI

    filenameam=$(echo "$file" | sed 's/_RG$//')
    singularity exec /cvmfs/singularity.galaxyproject.org/s/r/sra-tools:3.0.3--h87f3376_0 fastq-dump --split-3 -O $TMPWD $filenameam

    singularity exec tabix_1_11_hdfd78af_0.sif bgzip -f ${filenameam}_1.fastq > ${file}_1.fastq
    singularity exec tabix_1_11_hdfd78af_0.sif bgzip -f ${filenameam}_2.fastq > ${file}_2.fastq
  else
      ##dl file
      report=$WD/filereport_read_run_PRJEB38061_tsv.txt

      result=$(grep "$file" "$report" | awk '{print $10}')
      result=${result%;*}
      result=${result%_*}

      wget -c ${result}_R1.fastq.gz -P $TMPWD/${file}_1.fastq.gz
      wget -c ${result}_R2.fastq.gz -P $TMPWD/${file}_2.fastq.gz

      # Specify the folder name
      folder=$TMPWD/${file}_1.fastq.gz

      # Get the file name from the folder
      file2=$(ls "$folder")

      # Rename the file
      mv "$folder/$file2" "$file2"

      # Remove the folder
      rm -r "$folder"

      ##rename file again
      mv "$file2" "$folder"

    ## same for read 2
      folder=$TMPWD/${file}_2.fastq.gz

      file2=$(ls "$folder")

      mv "$folder/$file2" "$file2"

      rm -r "$folder"

      mv "$file2" "$folder"
  fi
    
    ##make the new file with the fasta file name

    ls -d $PWD/* |grep gz > lista.txt

    singularity exec kmc_3_2_1_h95f258a_1.sif kmc -k55 -ci1 -fq @lista.txt $file $TMPWD

    ## copy back some files

    rm *.fastq.*
    ## make bloom filter
    echo "Making bloom filter" 
    singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyperTools makeBloom -k ${file}

done


for value in "${selected_values[@]}"; do
    grep -w "$value" "$WD/bayertyper/samples.tsv" >> "tempsamples.tsv"
done

### get some files
cp $WD/bayertyper/combined_files.vcf.gz ./
cp $WD/data/Simon_Final2021_CHR.fasta ./
cp $WD/containers/bcftools_1_9_ha228f0b_4.sif ./
cp $WD/chr_ploidy.txt ./

##identify variant cluster 
singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyper cluster -v combined_files.vcf.gz -s tempsamples.tsv -g Simon_Final2021_CHR.fasta -o cluster_file

##genotype 

##count the number of cluster file
CLUSTNB=$(ls cluster_file_unit* | wc -l)

if [ "$CLUSTNB" -eq 1 ]; then
    # Perform operation for value = 1
    echo "1 cluster file"

    singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyper genotype -v cluster_file_unit_1/variant_clusters.bin -c cluster_file_cluster_data -s tempsamples.tsv --chromosome-ploidy-file chr_ploidy.txt -g Simon_Final2021_CHR.fasta -o genotyped_file_${SLURM_ARRAY_TASK_ID}

    ## save results
    cp genotyped_file* $WD/bayertyper/

elif [ "$CLUSTNB" -eq 2 ]; then
    # Perform operation for value = 2
    echo "2 cluster file"

    singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyper genotype -v cluster_file_unit_1/variant_clusters.bin -c cluster_file_cluster_data -s tempsamples.tsv --chromosome-ploidy-file chr_ploidy.txt -g Simon_Final2021_CHR.fasta -o genotyped_file_${SLURM_ARRAY_TASK_ID}_1
    singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyper genotype -v cluster_file_unit_2/variant_clusters.bin -c cluster_file_cluster_data -s tempsamples.tsv --chromosome-ploidy-file chr_ploidy.txt -g Simon_Final2021_CHR.fasta -o genotyped_file_${SLURM_ARRAY_TASK_ID}_2

    ##concat results
    bcftools_1_9_ha228f0b_4.sif bcftools concat -O z -o genotyped_file_${SLURM_ARRAY_TASK_ID}.vcf.gz genotyped_file_${SLURM_ARRAY_TASK_ID}_1.vcf genotyped_file_${SLURM_ARRAY_TASK_ID}_2_vcf

    ## save results
    cp genotyped_file*.vcf.gz $WD/bayertyper/

elif [ "$CLUSTNB" -eq 3 ]; then
    # Perform operation for value = 3
    echo "3 cluster file"

    singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyper genotype -v cluster_file_unit_1/variant_clusters.bin -c cluster_file_cluster_data -s tempsamples.tsv --chromosome-ploidy-file chr_ploidy.txt -g Simon_Final2021_CHR.fasta -o genotyped_file_${SLURM_ARRAY_TASK_ID}_1
    singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyper genotype -v cluster_file_unit_2/variant_clusters.bin -c cluster_file_cluster_data -s tempsamples.tsv --chromosome-ploidy-file chr_ploidy.txt -g Simon_Final2021_CHR.fasta -o genotyped_file_${SLURM_ARRAY_TASK_ID}_2
    singularity exec bayestyper_1_5_h176a8bc_0.sif bayesTyper genotype -v cluster_file_unit_3/variant_clusters.bin -c cluster_file_cluster_data -s tempsamples.tsv --chromosome-ploidy-file chr_ploidy.txt -g Simon_Final2021_CHR.fasta -o genotyped_file_${SLURM_ARRAY_TASK_ID}_3

    ##concat results
    bcftools_1_9_ha228f0b_4.sif bcftools concat -O z -o genotyped_file_${SLURM_ARRAY_TASK_ID}.vcf.gz genotyped_file_${SLURM_ARRAY_TASK_ID}_1.vcf genotyped_file_${SLURM_ARRAY_TASK_ID}_2_vcf genotyped_file_${SLURM_ARRAY_TASK_ID}_3_vcf

    ## save results
    cp genotyped_file*.vcf.gz $WD/bayertyper/
else
    ### If more cluster, don't run
    echo "More than 3 cluster file, terminating"
fi


fi ##end of step2



