args <- commandArgs(trailingOnly = TRUE)

library(SeqArray)
library(survival)
library(parallel)

load('~/hap/data/ukb.bca.pheno.rdata')
source("~/hap/scripts/gdsCoxFit4.R")

gdsfn  <- args[1]
w      <- as.numeric(args[2])
s      <- as.numeric(args[3])
type   <- args[4]
fcut1   <- as.numeric(args[5])
fcut2   <- as.numeric(args[6])
ncore  <- as.numeric(args[7])
outdir <- args[8]
bsize  <- as.numeric(args[9])

gds <- seqOpen(gdsfn)
chrid <- as.numeric(seqGetData(gds, 'chromosome')[1])
vid <- seqGetData(gds, 'variant.id')
set.list <- splitSet(vid, w=w, s=s, type=type)
seqClose(gds)

logfile <- paste(outdir,'/chr',chrid,'.log',sep='')
write.table(paste('chromosome ', chrid, ' has ', length(set.list), ' blocks...', sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)

system(paste('mkdir ', outdir, '/chr', chrid, sep=''))
blocks <- split(1:length(set.list), ceiling(seq_along(set.list)/bsize))
write.table(paste('scan chromosome ', chrid, ' in ', length(blocks), ' segments...', sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
write.table(paste('each segment has ', bsize, ' blocks', sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)

for(j in 1:length(blocks)){
    bset <- set.list[blocks[[j]]]
    write.table(paste('Chromosome ', chrid, ' >>>> ', 'segment ', j, ' calcuate MAFs in BCa cases', sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
    fit <- mclapply(bset, gdsCoxFit, gdsfn=gdsfn, pheno=ukb.bca.pheno, 
		       yy=c("before.analysis.end", "before.analysis.diag"), xx=paste('pc',1:10,sep=''), type=type, 
               frq.min=fcut1, frq.max=fcut2, maf_sids=ukb.bca.pheno[ukb.bca.pheno$before.analysis.diag==1, ]$f.eid,
               mc.cores = ncore)
    fit <- do.call('rbind', fit)
    save(fit, file=paste(outdir,'/chr', chrid, '/b', j, '.rdata',sep=''))
    write.table(paste('chromosome ', chrid, ' >>>> ', 'segment ', j, ' completed.', sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
    rm(fit)
}

## collect results
write.table('collecting results >>>>>>>', file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
res <- list()
for(j in 1:length(blocks)){
    write.table(paste('segment ', j, sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
    file.tmp <- paste(outdir,'/chr', chrid, '/b', j, '.rdata', sep='')
    if(file.exists(file.tmp)){
          load(file.tmp)
          res[[j]] <- fit
          rm(fit)
    }else{
          write.table(paste(file.tmp, ' does not exist..', sep=''), file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)
    }
}
res <- do.call('rbind',res)
save(res, file=paste(outdir,'/', type, '_w', w, '_s', s, '_chr', chrid, 'mafCaseOnly.rdata', sep=''))

## clean up
system(paste('rm -r ', outdir, '/chr', chrid, sep=''))

write.table('analysis completed...', file=logfile, sep='\n', col.names=F, row.names=F, quote=F, append=T)