## General information

This Directory contains script used in [paper name]



the subfolder **bayestyper_SV_calling**; **SV_PCA_calling** and **slim** are used to generate data, see their own readme in the subfolder for more information.
The subfolder **Others** contains others scripts for analysis and plotting, with content detailled below.



## General usage



The script <ins>get_vcf_stat.sh</ins> is a tool used in many others scripts, including those used to generate data. It assemble some useful commands to extract information from vcf files.





## Analysis of Bayestyper_SVs



- <ins>Sv_bayestyper_metrics.R</ins> is used to analyse basic stats of bayestyper_SV, such as their size, frequency, overlap With gene annotations (Figure 1, Supplementary figure 2)


- <ins>repeatSVanalysis.sh</ins> is used to analyse SV content around SVs. Some R scripts codes called within this script can be found in <ins>SV_TE_overlap.R</ins>, which also contains the code to analyse and plot the TE content around SV (Figure 2)


- <ins>ibs_haplotype_bayestyperSV.sh</ins> is used to run ibs analysis on bayestyper_SVs, and plots can be made using <ins>plot_bayestyper_ibs.R</ins>.


- <ins>Pca_outlier.R</ins> is used to make a PCA using besytyper_SVs which helped to identify the outlier sample (Supplementary figure 3).



## Analysis of lostruct_SVs


Most plots of Figure 3 and 4 are done using scripts in the **SV_PCA_calling** directory.



- <ins>get_sample_read_depth.sh</ins> is used to get reads depth for each individual, and <ins>read_depth.R</ins> is used to format into a Matrix. This Matrix is then used in scripts in the **SV_PCA_calling** directory.


- <ins>frequency_SV_plot.R</ins> is used to plot frequency of lostruct_SVs (and bayestyper_SVs) (Figure 5).


- <ins>xpehh_overlap.R</ins> is used to investigate overlap of lostructs_SVs and bayestyper_SVs with xpehh data from *Buso, P., Diblasi, C., Manousi, D., Kwak, J. S., Vera-Ponce de Leon, A., StenlC8kk, K., ... & Saitou, M. (2025). Parallel Selection in Domesticated Atlantic Salmon from Divergent Founders Including on Whole-Genome Duplication-derived Homeologous Regions. Genome Biology and Evolution, 17(4), evaf063.* (Supplementary figure 6)


- <ins>ibs_haplotype_LostructSV.sh</ins> is used to run ibs analysis on lostruct_SVs, and plots are made using <ins>plot_lostruct_ibs.R</ins> (Figure 6 B,C and Supplementary figure 8).


- <ins>TE_recurrent_Svs.R</ins> is used to analyze the TE content around reccurent and non recurrent lostruct_SVs (Supplementary figure 9).



## Data

The **Others/Data** folder include the main data generated and used in others scripts:

- SV_info_table_farmed_Europe.txt and Bayestyper_Sv_info_table_farmed_europe.txt: respectively manta and bayestyper SV list and coordinates

- df_nb_repeat_global.txt and df_nb_repeat_global_extended.txt: the number of TEs of different types repectively within SVs or surrounding SVs.

- overlap_SV_annot_v2.txt: genome features overlapped by bayestyper_SVs.

- interesting_region_bis.txt: the list of the 25 lostructs SV regions and their coordinate.






