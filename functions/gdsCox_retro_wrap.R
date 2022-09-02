#************************************************
# step up
#************************************************
library(SeqArray)
library(survival)
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
## phenotype file
phenfile <- as.character(args[1])
## haplotyep windows 
setfile  <- as.character(args[2])
## Genotype GDS file 
gdsdir   <- as.character(args[3])
## output direcotry
outdir   <- as.character(args[4])
## number of CPU cores
ncore    <- as.numeric(args[5])

## Cox regression function
source("gdsCoxFit6.R")

#************************************************
# Loading data
#************************************************
load(phenfile)
load(setfile)
fname <- gsub('.rdata','', basename(setfile))
chrid <- as.numeric(substr(strsplit(fname,'[.]')[[1]][1],4,5))
gdsfn <- gsub('xxx', paste0('chr', chrid), gdsdir)
## creat temp directory
system(paste0('mkdir ', outdir, '/', fname))
## creat log file
logfile <- paste0(outdir,'/', fname,'.log')
write.table(paste('chromosome ', chrid, ': ', length(bsets), ' chunks...', sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)

#************************************************
# Testing Haplotypes by Cox regression
#************************************************
## testing each window blocks
for(j in 1:length(bsets)){
    bset <- bsets[[j]]
    fit <- mclapply(bset, gdsCoxFit, gdsfn=gdsfn, pheno=pheno, yy=c("before.analysis.end", "before.analysis.diag"), 
                    xx=paste('pc',1:10,sep=''), type='hap', frq.min=0.001, frq.max=0.500, 
                    maf_sids=pheno[pheno$before.analysis.diag==1, ]$f.eid, mc.cores = ncore)
    fit <- do.call('rbind', fit)
    save(fit, file=paste0(outdir,'/',fname, '/b', j, '.rdata'))
    write.table(paste0('chromosome ', chrid, ' >>>> ', 'chunk ', j, ' completed.'), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
    rm(fit)
    print(j)
}

## collect results
write.table('collecting results >>>>>>>', file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
res <- list()
for(j in 1:length(bsets)){   
    write.table(paste0('reading chunk ', j), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
    file.tmp <- paste0(outdir, '/', fname, '/b', j, '.rdata')
    if(file.exists(file.tmp)){
          load(file.tmp)
          res[[j]] <- fit
          rm(fit)
    }else{
          write.table(paste0(file.tmp, ' does not exist..'), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
    }
}
res <- do.call('rbind',res)
# save data
save(res, file=paste0(outdir, '/', fname, '.res.rdata'))

## clean up
system(paste0('rm -r ', outdir, '/', fname))
write.table('analysis completed...', file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
