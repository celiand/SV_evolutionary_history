##General information

These scripts are used to generate the lostruct_SV set. 


##Content

Outliers region where identified using the lostruct packages and the built in script, see https://github.com/petrelharp/local\_pca 



Detected outlier regions where then recorded in a files and a grouping was done using kmeans (see Methods); to generate a file With all regions and assign each individual to one group.



- <ins>fst_region_bycluster.sh</ins> is then used to plot FST of all possible contrast (between all Groups) within an outlier region. <ins>fst_cluster_process.R</ins> and <ins>plot_FST_bycluster.R</ins> are used with this script. 





Then, all plots are manually investigated and regions containing high FST blocks are recorded, with an approximate boundary. 



- <ins>Sliding_pca_region_general_investigation_narrow_version.sh</ins> is then used to compute and plots differents metrics, by taking as input the chromosome name, start and end of the region of interest. This can also be used to refine boundary location using the haplotype plot generated.

