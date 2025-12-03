## General information

This Directory contains script used in [paper name]



the subfolder **bayestyper_SV_calling**; **SV_PCA_calling** and **slim** are used to generate data, see their own readme in the subfolder for more information.
The subfolder **Others** contains others scripts for analysis and plotting, with content detailled below.



## general usage



The script <ins>get_vcf_stat.sh</ins> is a tool used in many others scripts, including those used to generate data. It assemble some useful commands to extract information from vcf files.





## analysis of Bayestyper_SVs



- <ins>Sv_bayestyper_metrics.R</ins> is used to analyse basic stats of bayestyper_SV, such as their size, frequency, overlap With gene annotations (Figure 1)



- <ins>repeatSVanalysis.sh</ins> is used to analyse SV content around SVs. Some R scripts codes called within this script can be found in <ins>SV_TE_overlap.R</ins>, which also contains the code to analyse and plot the TE content around SV (Figure 2)



## analysis of lostruct_SVs



Most plots of Figure 3 and 4 are done using scripts in the SV\_PCA\_calling drectory



\- get\_sample\_read\_depth.sh is used to get reads Depth for each individual, and read\_depth.R is used to format into a Matrix. This Matrix is then used in scripts in the SV\_PCA\_calling drectory



\- frequency\_SV\_plot.R is used to plot frequency of PCA\_SV (and bayestyper\_SV) (Figure 5)



