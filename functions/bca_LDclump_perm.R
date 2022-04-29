args <- commandArgs(trailingOnly = TRUE)

library(parallel)
library(SeqArray)

source('~/hap/scripts/gdsGeno.R')
chrid <- as.numeric(args[1])
ncore <- as.numeric(args[2])
bsize <- as.numeric(args[3])
fn.phen <- args[4]
fn.res  <- args[5]
fn.out  <- args[6]
gdsdir <- '~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr'

load(fn.phen)
gds <- seqOpen(paste(gdsdir, 1, '.gds',sep=''), allow.duplicate=TRUE)
seqSetFilter(gds, sample.id=pheno$f.eid)
sid <- seqGetData(gds, 'sample.id')
pheno <- pheno[match(sid, pheno$f.eid), ]
table(pheno$f.eid==sid)
seqClose(gds)

geno_extract <- function(idx, res, gdsdir, sid){
      tt <- as.numeric(strsplit(res[idx, ]$block,'[_]')[[1]])
      gdsfn <- paste(gdsdir, tt[1], '.gds',sep='')
      gg <- gdsGeno(tt[2]:tt[3], gdsfn, sid, type='hap', as.character(res[idx,]$hap))
      return(gg)
}

get_ped <- function(idx, geno){
       tt <- sapply(as.character(geno[,idx]), function(x){switch(x, 'NA'='0 0', '0'='1 1', '1'='1 2', '2'='2 2')})
       print(idx)
       return(tt)
}

load(fn.res)
assign('res', get(gsub('.rdata','', basename(fn.res))))
top <- res[res$chr==chrid, ]

if(nrow(top)>0){
   sets <- split(1:nrow(top), ceiling(seq_along(1:nrow(top))/bsize))
   for(i in 1:length(sets)){
      top.i <- top[sets[[i]], ]
      gg <- mclapply(1:nrow(top.i), geno_extract, res=top.i, gdsdir=gdsdir, sid=pheno$f.eid, mc.cores=ncore)
      gg <- do.call('cbind', gg)
      rownames(gg) <- pheno$f.eid
      colnames(gg) <- top.i$block

      ped <- mclapply(1:ncol(gg), get_ped, geno=gg, mc.cores=ncore)
      ped <- do.call('cbind', ped)
      nn <- nrow(ped)
      tmp <- data.frame(fid=rep(0,nn), iid=pheno$f.eid, pid=0, mid=0, sex=pheno$sex, phen=-9)
      ped <- data.frame(tmp, ped)
      map <- data.frame(chr=as.numeric(top.i$chr), snp=top.i$hid, cm=0, pos=top.i$start)
      options(scipen=999)
      write.table(ped, file=paste(fn.out, chrid,'.set',i,'.ped',sep=''), sep='\t',col.names=F, row.names=F, quote=F)
      write.table(map, file=paste(fn.out, chrid,'.set',i,'.map',sep=''), sep='\t',col.names=F, row.names=F, quote=F)
      print(i)
      }
}
