## **Usage** ##

### 1. Finding the SNVs with wide posterior distribution: ###
You can use CI_filter.py for this purpose. This tool calculates the quantile based credible interval for the prevalence of each SNV. The arguments are as follows:

|Argument| Description|
|----|-----------|
| -i | size of the quantile-based credible interval, CI, for the prevalence of each SNV. i must be a number between 0 and 100. (see below for more details on how to obtain a quantile-based credible interval) |
| -r | the fraction of samples allowed to have CIs that exceed a maximum designated length. CI_filter.py considers designated lengths between 0.1 and 0.9 with 0.1 increments. |
| -c | path to the PyClone config file. Please note that if you have copied PyClone results from their original location, you need to update working_dir and trace_dir parameters in the config file. | 
| -o | path to the analysis results. | 
| -b | burnin for the posterior distribution. |  
  
  
#### **Output files:** ####
CI_filter.py generates 10 files in the path that you are specified in -o argument.

1. summary_*.csv contains the fraction of SNVs with i% CI spanning longer than each maximum designated length in more than (100 x r)% of the samples.
2. For each maximum length X there is a file named snps_list_CI_X_node_rej_r.csv. This file contains a comma separated list of the SNVs with CIs that exceed the maximum allowed CI length in more than (100 x r)% of the samples.

### **2. Removing the SNVs with non-informative posterior distributions (for short non-informative SNVs)** ###
After finding the SNVs with large credible intervals in the previous stage using CI_filter.py, the user can apply pyclone_files_updater.py to remove those SNVs from all PyClone files and re-cluster the remaining SNVs.

User needs to decide the set of SNVs to be removed by analyzing the fraction of removable SNVs and the corresponding maximum allowable length of the CIs and use the corresponding file with the list of removable SNVs from stage 1.

The arguments for pyclone_files_updater.py are as follows:

| Argument| Description| 
|----------|-------------|
| -c | path to the PyClone config file. Please note that if you have copied PyClone results from their original location, you need to update working_dir and trace_dir parameters in the config file. |
| -s | path to the csv file containing the list of the removable SNVs. This file is generated in stage 1. |
| -f | path to the new cluster file. pyclone_files_updater.py clusters the remaining SNVs after removing the non-informative SNVs. |
| -b | burnin for the posterior distribution. |  

Please note that after removing the non-informative SNVs, pyclone_file_updater.py generates new files and keeps the original files by adding .originalto the end of the original files before removing those SNVs.

### **3. Significance analysis of change of prevalence of clusters in multi-sample experiments** ###
In multi-sample experiments, prevalence of some clusters may change between samples. It is important to know how significant those changes are. To this end, we can measure for each cluster the overlap of the credible intervals between every pair of samples. If the CIs don’t overlap it means the change in prevalence is statistically significant.

interval_analyser.py can be used for this purpose. The arguments are as follows:

| Argument| Description| 
|----------|-------------|
| -c | path to the PyClone config file. |
| -f | path to the PyClone cluster file. You need to use the new cluster file, generated in stage 1. |
| -o | path to the analysis results. |
| -b | burnin for the posterior distribution. |
| -i | size of the quantile-based credible interval, CI, for the prevalence of each SNV. i must be a number between 0 and 100. |

Output file contains 8 columns: cluster id, x1 and x2 (the start and end of the credible interval of the first sample) y1 and y2 (the start and end of the credible interval in the second sample), the mean of the posterior samples of prevalence in the first and second sample and result of the credible interval overlap test. In the last column 0 means no overlap and 1 means overlap of the credible intervals.