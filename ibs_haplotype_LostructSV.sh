#!/bin/bash

#SBATCH --cpus-per-task=4            # 1 core(CPU)
#SBATCH --job-name=haplotypeibs  # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. 
#SBATCH --output=haplotypeibs.txt  #output

WD=Data/

DD=Datadir/



eval "$(micromamba shell hook --shell bash)"

micromamba activate $HOME/toolsenv




## loop over all Sv
tail -n +2 "Data/interesting_region_bis.txt" | sed 's/"//g' | while read -r chrom oldcode start end newcode length; do
    
module load PLINK/2.00a2.3_x86_64

# paste the values
region="${chrom}_${start}_${end}"

#bcftools view -r ${chrom}:${start}-${end} -O z -o $DD/sprgi_vcffiles/PCA_${chrom}_${start}_${end}_int_phased.vcf.gz Data/subset368_SNP_phased.vcf.gz

#tabix -p vcf $DD/sprgi_vcffiles/PCA_${chrom}_${start}_${end}_int_phased.vcf.gz


# We use plink2 to read the VCF and generate PLINK binary files (.bed/.bim/.fam).
plink2 \
  --vcf $DD/sprgi_vcffiles/PCA_${chrom}_${start}_${end}_int_phased.vcf.gz \
  --double-id \
  --set-all-var-ids @:#:\$r:\$a \
  --max-alleles 2 \
  --allow-extra-chr \
  --make-bed \
  --new-id-max-allele-len 1000 \
  --out $WD/ibsfolder/${chrom}_${start}_${end}_rg

module load PLINK/1.9b_6.17-x86_64

  # Here we use the “classic” plink (v1.9).
# --distance square 1-ibs produces a square matrix of pairwise distances.
# “1-IBS” means the proportion of sites that are *not identical-by-state*.
# The output is:
#   region_ibs.mdist     : distance matrix
#   region_ibs.mdist.id  : sample IDs (order of rows/columns)
plink \
  --bfile $WD/ibsfolder/${chrom}_${start}_${end}_rg \
  --distance square 1-ibs \
  --allow-extra-chr \
  --out $WD/ibsfolder/${chrom}_${start}_${end}_rg_ibs

done