qnormRT.R
==========
An R script to quantile-normalize Repli-chip and Repli-seq data sets.

####Usage
qnormRT( RepliSeqBedGraphs , RepliChipFiles, cores="max" )

####Arguments
RepliSeqBedGraphs - a character vector of length 1 or more containing file paths/names of Repli-seq data in bedGraph format.

RepliChipFiles - a character vector of length 1 or more containing file paths/names of Repli-chip data with the first column containing chromosome names, the second column containing chromosome coordinates, and additional columns containing replication timing scores, with each column corresponding to a sample.

cores - the number of threads to utilize. Default is "max" which sets "cores" to 1 less than the total available cores.

###Details
qnormRT uses Repli-seq and Repli-chip data to generate a pool of replication timing scores which is subsequently used to sample scores for quantile normalization of the input replication timing profiles. The approach of pooling and subsampling scores allows for the quantile-normalization of data sets with different numbers of data points. The output of qnormRT is a new bedGraph of quantile-normalized data for each input sample. This function expects that both Repli-chip and Repli-seq profiles have similar, zero-centered distributions with similar IQRs.

####Value
qnormRT creates a bedGraph file in the working directory for each sample corresponding to a quantile-normalized replication timing profile, and returns a character vector of file names for each bedGraph created.