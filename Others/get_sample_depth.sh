#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --output=logs/mosdepth_%A_%a.out
#SBATCH --error=logs/mosdepth_%A_%a.err
#SBATCH --array=0-35%13
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --partition=hugemem-avx2

# Conda 
eval "$(conda shell.bash hook)"
conda activate ~/.conda/envs/myenv

set -euo pipefail

# list all samples
SAMPLE_LIST=(
  MA1_RG.bam
  MA21_RG.bam
  MA22_RG.bam
  MA26_RG.bam
  MA28_RG.bam
  MA31_RG.bam
  MA5_RG.bam
  MA6_RG.bam
  MA8_RG.bam
  MA9_RG.bam
  PC10_RG.bam
  PC11_RG.bam
  PC12_RG.bam
  PC14_RG.bam
  PC19_RG.bam
  PC21_RG.bam
  PC24_RG.bam
  PC26_RG.bam
  PC7_RG.bam
  PC8_RG.bam
  SP16_RG.bam
  SP21_RG.bam
  SP22_RG.bam
  SP23_RG.bam
  SP27_RG.bam
  SP28_RG.bam
  SP29_RG.bam
  SP30_RG.bam
  SP31_RG.bam
  SP8_RG.bam
  VF19_RG.bam
  VF22_RG.bam
  VF26_RG.bam
  VF28_RG.bam
  VF29_RG.bam
  VF3_RG.bam
  VF7_RG.bam
  VF8_RG.bam
)

BAM="${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}"
BAM_PATH="mydirectory/$BAM"
SAMPLE=$(basename "$BAM_PATH" .bam)

OUTDIR="Data/mosdepth_results/${SAMPLE}"
mkdir -p "$OUTDIR"


# ==== 1. CHROM_LIST ====
CHROM_LIST=($(seq -f "NC_059%03g.1" 442 470))  
SSA_MAP=($(seq -f "ssa%02g" 1 29))  

# ==== 2. mosdepth ====
for i in "${!CHROM_LIST[@]}"; do
  CHROM="${CHROM_LIST[$i]}"
  SSA_NAME="${SSA_MAP[$i]}"
  PREFIX="${OUTDIR}/${SSA_NAME}"

  echo "[$SAMPLE] Processing $CHROM → $SSA_NAME"

  mosdepth \
    --no-per-base \
    --fast-mode \
    --flag 1796 \
    --by 5000 \
    --threads ${SLURM_CPUS_PER_TASK} \
    --chrom ${CHROM} \
    "${PREFIX}" \
    "${BAM_PATH}"
done
