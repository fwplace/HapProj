gds_fit <- function(set.idx, gdsfn, pheno, yy, xx, type, frq.min, frq.max, method, get.lrt=FALSE, mod0.stat=NA, genetic.model='add', maf_sids, haps){ 

## match phenotype and genotype data
gds <- seqOpen(gdsfn, allow.duplicate=TRUE)
sid <- intersect(seqGetData(gds, 'sample.id'), pheno$sid)
seqSetFilter(gds, sample.id=sid, variant.id=set.idx)
pheno <- pheno[match(sid, pheno$sid), ]
vid <- seqGetData(gds, 'variant.id')

## Testing SNPs
if(type=='snp'){
    geno <- seqGetData(gds, '$dosage_alt')
    if(sum(vid!=set.idx)>0){
       order.idx <- sapply(vid, function(x){which(set.idx==x)})
       geno <- geno[, order.idx]
    }
    if(missing(maf_sids)){
           eaf <- colMeans(geno, na.rm=T)/2
    }else{
           sid.idx <- which(pheno$sid%in%maf_sids)
           eaf <- colMeans(geno[sid.idx, ], na.rm=T)/2
    }
    maf <- eaf
    maf[which(eaf>0.5)] <- 1-eaf[which(eaf>0.5)]
    vidx <- which(maf>=frq.min&maf<frq.max)
	if(length(vidx)>0){    
        sid.idx <- which(pheno[, yy[length(yy)]]==1)
        if(length(vidx)==1){
            geno2 <- as.numeric(geno[, vidx])
            eaf.case <- mean(geno2[sid.idx],na.rm=T)/2
            eaf.conc <- mean(geno2[-sid.idx],na.rm=T)/2
            eaf <- mean(geno2, na.rm=T)/2
            fit <- get_pva(z=geno2, nz=1, y=pheno[, yy], x=pheno[, xx], method=method, get.lrt=get.lrt, mod0.stat=mod0.stat, genetic.model=genetic.model)
        }else{
            geno2 <- lapply(vidx, function(x){as.numeric(geno[, x])})
            eaf.case <- as.numeric(sapply(geno2, function(x){mean(x[sid.idx],na.rm=T)/2}))
            eaf.conc <- as.numeric(sapply(geno2, function(x){mean(x[-sid.idx],na.rm=T)/2}))
            eaf <- as.numeric(sapply(geno2, function(x){mean(x,na.rm=T)/2}))
            fit <- do.call('rbind', lapply(geno2, get_pva, nz=1, y=pheno[, yy], x=pheno[, xx], method=method, get.lrt=get.lrt, mod0.stat=mod0.stat, genetic.model=genetic.model))
        }
        res <- data.frame(vid=seqGetData(gds, 'variant.id')[vidx], ea=seqGetData(gds, '$alt')[vidx], nea=seqGetData(gds, '$ref')[vidx], eaf.case=eaf.case, eaf.conc=eaf.conc, eaf=eaf, fit)
        rownames(res) <- res$vid
	}else{
	      print('No SNP was tested, please check your data!')
          res <- NULL
	}
seqClose(gds)
}

## Testing Haplotypes 
if(type=='hap'){
    hh <- seqGetData(gds, 'genotype')
    h1 <- hh[1,,]
    h2 <- hh[2,,]
    ## if vid does not match set.idx
    if(sum(vid!=set.idx)>0){
       order.idx <- sapply(vid, function(x){which(set.idx==x)})
       h1 <- h1[, order.idx]
       h2 <- h2[, order.idx]
    }
    h1 <- apply(h1, 1, paste, collapse="")
    h2 <- apply(h2, 1, paste, collapse="")
    h1 <- paste('h',h1,sep='')
    h2 <- paste('h',h2,sep='')
    hh <- c(h1, h2)
    if(missing(haps)) haps <- unique(hh)
    nn <- length(hh)
    if(missing(maf_sids)){
          eaf <- as.numeric(sapply(haps, function(x){sum(hh==x)/nn}))
    }else{
          sid.idx <- which(pheno$sid%in%maf_sids)
          hh2 <- hh[c(sid.idx, nrow(pheno)+sid.idx)]
          nn2 <- length(hh2)
          eaf <- as.numeric(sapply(haps, function(x){sum(hh2==x)/nn2}))
    }
    maf <- eaf
    maf[which(eaf>0.5)] <- 1-eaf[which(eaf>0.5)]
    hidx <- which(maf>=frq.min&maf<frq.max)
    if(length(hidx)>0){
       sid.idx <- which(pheno[, yy[length(yy)]]==1)
       hh2 <- hh[c(sid.idx, nrow(pheno)+sid.idx)]
       nn2 <- length(hh2)
       eaf.case <- as.numeric(sapply(as.vector(unlist(haps[hidx])), function(x){sum(hh2==x)/nn2}))
       hh2 <- hh[-c(sid.idx, nrow(pheno)+sid.idx)]
       nn2 <- length(hh2)
       eaf.conc <- as.numeric(sapply(as.vector(unlist(haps[hidx])), function(x){sum(hh2==x)/nn2}))
       eaf <- as.numeric(sapply(as.vector(unlist(haps[hidx])), function(x){sum(hh==x)/nn}))
       hinfo <- data.frame(hap=haps[hidx], eaf.case=eaf.case, eaf.conc=eaf.conc, eaf=eaf)
       hinfo <- hinfo[order(-hinfo$eaf), ]
       hinfo <- data.frame(hid=paste('h',1:nrow(hinfo), sep=''), hinfo)
       chr <- as.numeric(seqGetData(gds, 'chromosome'))
       pos <- as.numeric(seqGetData(gds, 'position'))
       if(nrow(hinfo)==1){
              hmat <- as.numeric(h1==hinfo$hap)+as.numeric(h2==hinfo$hap)
              fit <- get_pva(z=hmat, nz=1, y=pheno[, yy], x=pheno[, xx], method=method, get.lrt=get.lrt, mod0.stat=mod0.stat, genetic.model=genetic.model)
       }else{
            hmat <- do.call('cbind', lapply(hinfo$hap, function(x){as.numeric(h1==x)})) + do.call('cbind', lapply(hinfo$hap, function(x){as.numeric(h2==x)}))
	        colnames(hmat) <- hinfo$hid
            hmat2 <- lapply(1:ncol(hmat), function(x){as.numeric(hmat[, x])})
            fit <- do.call('rbind', lapply(hmat2, get_pva, nz=1, y=pheno[, yy], x=pheno[, xx], method=method, get.lrt=get.lrt, mod0.stat=mod0.stat, genetic.model=genetic.model))
       }
       res <- data.frame(block=paste(chr[1], set.idx[1], set.idx[length(set.idx)], sep='_'), chr=chr[1], start=min(pos), end=max(pos), t0=length(haps), t1=nrow(hinfo), hinfo, fit)
    }else{
	       print('No Haplotype was tested, please check your data!')
           res <- NULL
    }
    seqClose(gds)
}
return(res)
}
