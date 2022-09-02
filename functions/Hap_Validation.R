args <- commandArgs(trailingOnly = TRUE)
library(SeqArray)
library(survival)
source('~/fwork/ProjHap/scripts/gds.calc.pva.R')
setobj <- as.character(args[1])
resfn  <- as.character(args[2])

gdsfn <- '~/fwork/ProjHap/data/dbgap28544/topmed/qced/dbgap.topmedr2.ukb.xxx.qced.hg19.gds'
pp <- read.table(file='~/fwork/ProjHap/data/dbgap28544/dbgap28544.pheno.txt', sep='\t', header=T)
load(setobj)

validation_test <- function(hdse, hcut.type, hcut.value, exempt, exempt.id, y, x){
fit <- list()
for(i in 1:length(hcut.value)){
    if(exempt==0){
          print('No exempted SNP....')
          print(paste0('testing ', hcut.type, ' of hdse <= ', hcut.value[i]))
          tt <- as.numeric(apply(hdse, 1, hcut.type)<=hcut.value[i])
    }
    if(exempt==1){
          print('1 exempted SNP....')
          print(paste0('testing ', hcut.type, ' of hdse <= ', hcut.value[i]))
          tt <- as.numeric(apply(hdse[, -exempt.id], 1, hcut.type)<=hcut.value[i]&hdse[, exempt.id]<1)
    }
    if(exempt==2){
          print('1 exempted SNP....')
          print(paste0('testing ', hcut.type, ' of hdse <= ', hcut.value[i]))
          t1 <- as.numeric(apply(hdse[, -exempt.id], 1, hcut.type)<=hcut.value[i]&apply(hdse[, exempt.id], 1, max)<1)
    }
    nn <- length(tt)/2
    z <- tt[1:(nn)] + tt[(nn+1):(2*nn)]
    fit1 <- fisher_p(z=z, y=y, mod='dom')
    fit2 <- logistic_p(z=z, nz=1, y=y, x=x, get.lrt=TRUE)
    fit[[i]]  <- data.frame(type=hcut.type, hcut=hcut.value[i], exempt=exempt, 
                            snpids=ifelse(exempt==0, '.', paste(exempt.id, collapse='/')),
                            frq0=1e4*mean(z[y==0])/2, frq1=1e4*mean(z[y==1])/2,
                            OR=fit2$or, Fisher=fit1$p, LRT=fit2$p.lrt)
}
fit <- do.call('rbind', fit)
return(fit)
}

resObj <- list()
for(i in 1:13){
  ## Extract DRIVE imputed data  
  snp <- obj[[i]]$snp[obj[[i]]$snp$validation==1, ]
  gds <- seqOpen(gsub('xxx',paste0('chr', snp$chr[1]), gdsfn), allow.duplicate=TRUE)
  sid <- intersect(seqGetData(gds, 'sample.id'), pp$sid)
  phen <- pp[match(sid, pp$sid), ]
  vv <- data.frame(vid=seqGetData(gds,'variant.id'), id=seqGetData(gds,'annotation/id'))
  vv$id <- gsub(':','.',vv$id)
  tmp <- intersect(vv$id, snp$id)
  vv <- vv[match(tmp, vv$id), ]
  seqSetFilter(gds, sample.id=sid, variant.id=vv$vid)
  hds  <- seqGetData(gds, "annotation/format/HDS")$data
  hds  <- rbind(hds[, seq(1,ncol(hds),2)], hds[, seq(2, ncol(hds), 2)])
  hdse <- abs(sweep(hds, 2, as.numeric(snp$allele)))
  seqClose(gds)

  ## validation using orginal calls
  r1 <- validation_test(hdse=hdse, hcut.type='max', hcut.value=0.5, 
                        exempt=0, exempt.id=NULL, y=phen[, 'BCa'], x=phen[, c('ageonset',paste0('pc',1:10))])

  ## extented analysis 1
  r2 <- validation_test(hdse=hdse, hcut.type='sum', hcut.value=c(seq(0, 0.01, 0.001),seq(0.02, 0.5,0.01)),
                       exempt=0, exempt.id=NULL, y=phen[, 'BCa'], x=phen[, c('ageonset',paste0('pc',1:10))])

  ## extended analysis 2
  r3 <- list()
  for(j in 1:ncol(hdse)){
      r3[[j]] <- validation_test(hdse=hdse, hcut.type='sum', hcut.value=c(seq(0, 0.01, 0.001),seq(0.02, 0.5,0.01)),
                       exempt=1, exempt.id=j, y=phen[, 'BCa'], x=phen[, c('ageonset',paste0('pc',1:10))])
  }

  ## Combined results 
  resObj[[i]] <- rbind(r1, r2, do.call('rbind', r3))
  rm(snp, sid, vv, phen, geno, hds, hdse, r1, r2, r3)
}
save(resObj, file=resfn)

