#!/bin/bash

#SBATCH --cpus-per-task=1            # 1 core(CPU)
#SBATCH --job-name=ibshaplotype  # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. 
#SBATCH --output=ibshaplotype.txt  #output


WD=Data/

resdir=$WD/splittedvcfbt

eval "$(micromamba shell hook --shell bash)"

micromamba activate $HOME/toolsenv

### use table that contains coordinates for all SVs
sed 's/"//g' "$WD/bayestyper_SV_info_table_Farmed_American.bed" | while read -r chrom start end id freq pop typesv; do
    
module load PLINK/2.00a2.3_x86_64

## make a vcf containing only the SV
bcftools view -r ${chrom}:${start}-${end} -O z -o $resdir/${chrom}_${start}_${end}_phased.vcf.gz /mnt/SCRATCH/cedi/phDSalmon/evolutionary_analysis/SV_evol/subset368_SNP_phased.vcf.gz

tabix -p vcf $resdir/${chrom}_${start}_${end}_phased.vcf.gz

# We use plink2 to read the VCF and generate PLINK binary files (.bed/.bim/.fam).
plink2 \
  --vcf $resdir/${chrom}_${start}_${end}_phased.vcf.gz \
  --double-id \
  --set-all-var-ids @:#:\$r:\$a \
  --max-alleles 2 \
  --allow-extra-chr \
  --make-bed \
  --new-id-max-allele-len 1000 \
  --out $WD/ibsfolder_bt/${chrom}_${start}_${end}_rg

module load PLINK/1.9b_6.17-x86_64

  # Here we use the “classic” plink (v1.9).
# --distance square 1-ibs produces a square matrix of pairwise distances.
# “1-IBS” means the proportion of sites that are *not identical-by-state*.
# The output is:
#   region_ibs.mdist     : distance matrix
#   region_ibs.mdist.id  : sample IDs (order of rows/columns)
plink \
  --bfile $WD/ibsfolder_bt/${chrom}_${start}_${end}_rg \
  --distance square 1-ibs \
  --allow-extra-chr \
  --out $WD/ibsfolder_bt/${chrom}_${start}_${end}_rg_ibs


done