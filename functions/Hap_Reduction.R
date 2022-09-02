args <- commandArgs(trailingOnly = TRUE)
library(SeqArray)
library(rpart)
library(survival)
source('~/fwork/ProjHap/scripts/gds.calc.pva.R')
setobj <- as.character(args[1])
r2cut1 <- as.numeric(args[2])
r2cut2 <- as.numeric(args[3])
msplit <- as.numeric(args[4])
resfn  <- as.character(args[5])
gdsfn <- '~/fwork/ProjHap/data/ukb.bca/gds/ukb.bca.hap.xxx.gds'
pp <- read.table(file='~/fwork/ProjHap/data/ukb.bca/gds/ukb.bca.pheno.txt', sep='\t', header=T)

load(setobj)
obj <- get(gsub('.rdata','', basename(setobj)))
for(i in 1:13){
    ## phased haplotype data
    snp <- obj[[i]]$snp
    gds <- seqOpen(gsub('xxx',paste0('chr', snp$chr[1]), gdsfn), allow.duplicate=TRUE)
    sid <- intersect(seqGetData(gds, 'sample.id'), pp$sid)
    seqSetFilter(gds, sample.id=sid, variant.id=snp$vid)
    phen <- pp[match(sid, pp$sid), ]
    geno <- seqGetData(gds, 'genotype')
    geno <- as.data.frame(rbind(geno[1,,], geno[2,,]))
    colnames(geno) <- snp$snp
    seqClose(gds)
    haps <- paste0('h', apply(geno, 1, paste, collapse=""))
    outcome <- as.numeric(haps==paste(c('h', snp$allele), collapse=''))
    
    ## Fitting rpart
    print('SNP Reduction')
    maf <- snp$comb.af
    maf[which(snp$comb.af>0.50)] <- 1-maf[which(snp$comb.af>0.50)]
    rsq <- as.numeric(snp$imp.rsq)
    rsq[is.na(rsq)] <- 0
    snp$qc <- as.numeric((maf>=0.01&rsq>=r2cut1)|(maf<0.01&rsq>=r2cut2))
    carrier_data <- data.frame(outcome=outcome, geno[, which(snp$qc==1)])
    rpart.fit=rpart(outcome ~., control = rpart.control(minsplit=msplit, cp=10^-10),data=carrier_data)
    snp$validation <-  as.numeric(snp$snp%in%rpart.fit$frame[rpart.fit$frame[,1]!="<leaf>",1])
    
    ## check reduced haplotype in the UKBB data
    print('Evaluate reduced SNPs in the UKBB data')
    vidx <- which(snp$validation==1)
    nn  <- nrow(geno)/2
    hh <- apply(geno[, vidx], 1, paste, collapse="")
    tt <- as.numeric(hh==paste0(snp$allele[vidx],collapse = ''))
    gg  <-  tt[1:nn] + tt[(nn+1):(2*nn)]
    fit <- cox_p(z=gg, nz=1, y=phen[, c('age.end', 'diag')], x=phen[, paste('pc',1:10,sep='')])
    tmp <- table(predict(rpart.fit),carrier_data$outcome)
    if(nrow(tmp)==2) tmp1 <- tmp[2, ]
    if(nrow(tmp)>2)  tmp1 <- colSums(tmp[-1, ])
    res <- data.frame(bid=snp$bid[1], hid=snp$hid[1], r2cut1=r2cut1, r2cut2=r2cut2, msplit=msplit,
                     snps=paste(c(nrow(snp), sum(snp$qc==1), sum(snp$validation==1)),collapse='|'),
                     risk0=as.numeric(tmp1)[1], risk1=as.numeric(tmp1)[2],
                     frq0=1e4*mean(gg[phen$diag==0])/2, 
                     frq1=1e4*mean(gg[phen$diag==1])/2, fit[, c(2,5)])
    
    obj[[i]]$snp <- snp
    obj[[i]]$ukb <- res
    print('DONE....')
}
save(obj, file=resfn)

