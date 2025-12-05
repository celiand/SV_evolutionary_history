#!/bin/bash
#SBATCH --cpus-per-task=24            # 1 core(CPU)
#SBATCH --job-name=Make_library   # sensible name for the job
#SBATCH --mem=100G                 # Default memory per CPU is 3GB.
#SBATCH --output=TE_library_f.txt  #output



WD=$HOME/TE_library

## command to run the script

## sbatch Make_TE_library.sh /mnt/project/Transpose/rewired/genomes/public/Ssal/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz Ssal


echo "** Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job. **"


GENOME=${1}
NAME=${2}

## check if genome is defined
if [ -z "$GENOME" ]; then
    echo "** No Genome provided. Terminating. **"
    exit
else
    echo "** Making TE library for genome $GENOME **"
fi

## check if name is defined
if [ -z "$NAME" ]; then
    echo "** No name provided. Terminating. **"
    exit
else
    echo "** Using name $NAME **"
fi

cd $WD

### First step: RepeatModeler ###
echo "** Building database **"
singularity exec $WD/containers/repeatmodeler_2_0_4_pl5321hdfd78af_0.sif BuildDatabase -name $NAME $GENOME

echo "** Making TE library with RepeatModeler **"
singularity exec $WD/containers/repeatmodeler_2_0_4_pl5321hdfd78af_0.sif RepeatModeler -database $NAME -LTRStruct -threads 24


### Then DeeptTE to classify unknown families ###

### Identification of Unknown families
echo "** identifying unknown TEs **"

FILE="$NAME-families.fa"

awk '/^>/{flag=0} /Unknown/{flag=1} flag' "$FILE" > unknown_TE_${NAME}.fasta

awk '/^>/{flag=0} /Unknown/{flag=1} !flag' "$FILE" > know_TE_${NAME}.fasta

module load Miniconda3 
eval "$(conda shell.bash hook)"
conda activate py36

dirname="resultsdeep$NAME"
mkdir "$dirname"

DeepTE/DeepTE.py -d $WD/tempdeep -o $WD/$dirname -i unknown_TE_${NAME}.fasta -sp M -m M

cat know_TE_${NAME}.fasta $WD/$dirname/opt_DeepTE.fasta > refined_${NAME}_families.fasta

### Then Cdhit to Remove redundancy ###

### Removing redundancy

singularity exec $WD/containers/cdhit_4_8_1_hdbcaa40_2.sif cd-hit-est -i refined_${NAME}_families.fasta -o refined_${NAME}_families_reduced.fasta -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500