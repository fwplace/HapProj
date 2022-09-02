args <- commandArgs(trailingOnly = TRUE)
library(SeqArray)
library(parallel)

gdsfn  <- as.character(args[1])
phenfn <- as.character(args[2])
setfn  <- as.character(args[3])
outfn <- as.character(args[4])
ncore  <- as.numeric(args[5])

gds_hap_count <- function(set.idx, gdsfn, pheno, yy){ 
## match phenotype and genotype data
gds <- seqOpen(gdsfn, allow.duplicate=TRUE)
sid <- intersect(seqGetData(gds, 'sample.id'), pheno$sid)
seqSetFilter(gds, sample.id=sid, variant.id=set.idx)
pheno <- pheno[match(sid, pheno$sid), ]
hh <- seqGetData(gds, 'genotype')
h1 <- hh[1,,]
h2 <- hh[2,,]
h1 <- apply(h1, 1, paste, collapse="")
h2 <- apply(h2, 1, paste, collapse="")
h1 <- paste('h',h1,sep='')
h2 <- paste('h',h2,sep='')
hh <- c(h1, h2)
haps <- as.vector(unlist(unique(hh)))
sid.idx <- which(pheno[, yy[length(yy)]]==1)
hh2 <- hh[c(sid.idx, nrow(pheno)+sid.idx)]
nn2 <- length(hh2)
eaf <- as.numeric(sapply(haps, function(x){sum(hh2==x)/nn2}))
chr <- seqGetData(gds, 'chromosome')[1]
pos <- seqGetData(gds, 'position')
frq <- c(0, 0.0001, 0.0005, 0.001, 0.01, 0.05)
res <- data.frame(block=paste0(chr, '_', set.idx[1], '_', set.idx[length(set.idx)]), 
                  size=length(set.idx), length=max(pos) - min(pos)+1, tot=length(haps), 
                  t1=sum(eaf>=frq[1]&eaf<frq[2]), t2=sum(eaf>=frq[2]&eaf<frq[3]),
                  t3=sum(eaf>=frq[3]&eaf<frq[4]), t4=sum(eaf>=frq[4]&eaf<frq[5]),
                  t5=sum(eaf>=frq[5]&eaf<frq[6]), t6=sum(eaf>=frq[6]))
seqClose(gds)
return(res)
}

load(setfn)
pheno <- read.table(file=phenfn, sep='\t', header=T)
res <- list()
for(j in 1:length(bsets)){
    bset <- bsets[[j]]
    fit <- mclapply(bset, gds_hap_count, gdsfn=gdsfn, pheno=pheno, yy='before.analysis.diag', mc.cores = ncore)
    res[[j]] <- do.call('rbind', fit)
    print(j)
    rm(fit)
}
res <- do.call('rbind', res)
write.table(res, file=outfn, sep='\t', col.names=T, row.names=F, quote=F)