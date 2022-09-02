# HapProj

This folder contains R functions and scripts used to conduct a genome-wide haplotype analysis using the UKBB GWAS data. 


1. Split the genome into smaller chunks given a fixed window size `w` (e.g., `w=20 variants`) and skip length `s` (e.g., `s=1 varaint`) 
 
```
Rscript ./functions/gds_splitSet.R [variant file] [w] [s] [type] [bsize] [jobsize] [outdir]
```   

* `variant file` is a R Dataframe containing two columns: chromosome and index of varaint in the genotype GDS file
* `w`, the number of consective varaints used to construct haplotypes 
* `s`, skip length of sliding windows
* `bsize`, the number of windows analyzed in a batch
* `jobsize`, the number of batches in a single job
* `outdir`, outout directory

2. 
