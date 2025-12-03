## General information 

These scripts are used to generate the bayestyper_SV set. 



## Content

- <ins>sv_detection_manta.sh</ins> is used to detect SV from individual bam file sample, and generate a vcf file containing detected SVs as output.



- <ins>filtering_and_merging_manta.sh</ins> is used to clean the obtained vcf and merge all samples together, to ensure that the same occurence of an SV are considered as a unique one accross samples.



- <ins>bayestyper_SV_genotyping.sh</ins> is used to genotype manta Sv in batch of 5 samples and generate a vcf for each batch.



- <ins>filtering_and_merging_bayestyper.sh</ins> is used to clean the vcf and merge them together to obtain an unique vcf.

