args <- commandArgs(trailingOnly = TRUE)

library(SeqArray)
library(survival)
library(parallel)

source('~/hap/scripts/gdsCoxFit6.R')

rep_wrap <- function(id, res, gdsdir, pheno, type){
   if(type=='snp'){
    fit <- gdsCoxFit(set.idx=as.numeric(res$vid[id]), 
                     gdsfn=paste(gdsdir, as.numeric(res$chr[id]), '.gds',sep=''), 
                     pheno=pheno,
                     yy=c('age_initial_assessment','age.end', 'diag'),
                     xx=paste('pc',1:10,sep=''),
                     type='snp', frq.min=0, frq.max=0.501)

   }

   if(type=='hap'){
     tt <- as.numeric(strsplit(res[id, ]$block,'[_]')[[1]])
     fit <- gdsCoxFit(set.idx=tt[2]:tt[3], 
                     gdsfn=paste(gdsdir, tt[1],'.gds',sep=''), 
                     pheno=pheno,
                     yy=c('age_initial_assessment','age.end', 'diag'),
                     xx=paste('pc',1:10,sep=''),
                     type='hap', frq.min=0, frq.max=0.501,
                     haps=as.character(res[id,]$hap))
   }
   return(fit)
}

chrid <- as.numeric(args[1])
ncore <- as.numeric(args[2])
fn.phen <- args[3]
fn.res  <- args[4]
fn.out  <- args[5]
gdsdir  <- '~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr'

load(fn.phen)
pheno <- pheno[pheno$before.analysis.diag==0, ]

load(fn.res)
assign('res', get(gsub('.rdata','', basename(fn.res))))
res <- res[res$chr==chrid, ]

if(nrow(res)>0){
   fit <- mclapply(1:nrow(res), rep_wrap, res=res, gdsdir=gdsdir, pheno=pheno, type='hap', mc.cores = ncore)
   fit <- do.call('rbind', fit)
   colnames(fit)[9:14] <- paste('pros.', colnames(fit)[9:14], sep='')
   res <- data.frame(res, fit[, c(9,11,12,14)])
   save(res, file=paste(fn.out, chrid, '.rdata',sep=''))
}
