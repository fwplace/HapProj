# HapProj

This folder contains R functions and scripts used to conduct a genome-wide haplotype analysis using the UKBB GWAS data. 


1. Split the genome into smaller chunks given a fixed window size `w` (e.g., `w=20 variants`) and skip length 
 

```
Rscript gds_splitSet.R [variant file] [w] [s] [type] [bsize] [jobsize] [oudir]
```   



2. Run haploty
