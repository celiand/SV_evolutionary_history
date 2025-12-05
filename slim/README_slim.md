## General information

These scripts are used to generate simulation data and analyse obtained results.



## Slim simulations

All SLiM scripts are located in the **SlimScript** directory. <ins>burn_in_pop</ins> generate a population of 4000 individuals after 10 000 generations, containing mutations, to be used as a strating point for further simulations.

The other scripts contains all simulated scenario, including directional selection (\*_selec) or heterozygote selection (\*_selec_het); as well as ancient origin or recurrent SVs.
All scripts generate a SV and then create 100 identical haplotypes, so the SV is less likely to be lost immediately. 


The scripts save a vcf containing 50 individuals, for each population.

The SLiM scripts can be called using <ins>sv_simu.sh</ins>.


The scripts are inspired by _Berdan, E. L., Blanckaert, A., Butlin, R. K., & Bank, C. (2021). Deleterious mutation accumulation and the long-term fate of chromosomal inversions. PLoS genetics, 17(3), e1009411._


## Analysis

After running script, <ins>Format_slim_result.sh</ins> to fromat vcf results in a single table, especially retrieving the type of each mutations (SNP or SV) and its frequency.

Further plots are done using <ins>analyze_svsimu.R</ins> (Figure 6A; Supplementary figure 7).


## Ibs analysis
For SVs generated using simulations, the Ibs analysis is done using <ins>ibs_haplotype_simudata.sh</ins> and <ins>run_ibsstat.R</ins>.

