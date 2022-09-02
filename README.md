# HapProj

This folder contains R functions and scripts used to conduct a genome-wide haplotype analysis using the UKBB GWAS data. 


1. Split the genome into smaller chunks given a fixed window size `w` (e.g., `w=20 variants`) and skip length `s` (e.g., `s=1 varaint`) 
 
```
Rscript functions/gds_splitSet.R [variant file] [w] [s] [type] [bsize] [jobsize] [outdir]
```   

* `variant file` is a R Dataframe containing two columns: chromosome and index of varaint in the genotype GDS file
* `w`, the number of consective varaints used to construct haplotypes 
* `s`, skip length of sliding windows
* `bsize`, the number of windows analyzed in a batch
* `jobsize`, the number of batches in a single job
* `outdir`, outout directory

2. Test all haplotypes for each small chunk using Cox regression

```
Rscript gdsCox_retro_wrap.R [phenfile] [setfile] [gdsdir] [outdir] [ncore]
```

* `phenfile` is a R Datafram contaning sample ids, disease status, age, and other covariates
* `setfile` is a R List returned in the step 1
* `gdsdir` is the genotype data saved in the [GDS format](https://bioconductor.org/packages/release/bioc/html/SeqArray.html)
* `outdir`, output directory
* `ncore`, the number of CPU cores

3. Merge all result files into a single file

```
Rscript gds.hap.sum.R [resPath] [outfn1] [outfn2]
```

* `resPath` is the output directory returned in the step 2
* `outfn1` outout all haplotypes analzyed
* `outfn2` output the most significant haplotype for each unique haplotype window
