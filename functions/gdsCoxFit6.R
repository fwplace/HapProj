# Fit Cox Regression Model
coxFit <- function(z, nz, y, x){
       
  if(nz==1){ 
         z <- as.numeric(z)
  }else{
         z <- as.matrix(z)
  }

  if(ncol(y)==2){
   	   colnames(y) <- c('t1', 'diag')
   	   dat <- cbind(z, y, x)
          fit <- tryCatch(coxph(Surv(t1, diag)~., data=dat), res <- rep(1, 5), error=function(e){res <- rep(NA, 5)})
   }else{
   	   colnames(y) <- c('t0', 't1', 'diag')
          dat <- cbind(z=z, y, x)
          fit <- tryCatch(coxph(Surv(t0, t1, diag)~., data=dat), res <- rep(1, 5), error=function(e){res <- rep(NA, 5)})
   }

   if(sum(is.na(res))==0){
   if(nz==1){
       	res <- summary(fit)$coef['z', ]
       	names(res) <- c('beta','hr', 'se', 'z', 'p')
      }else{
       	res <- summary(fit)$coef[colnames(z), ]
       	colnames(res) <- c('beta','hr', 'se', 'z', 'p')
      }
   }
   
   return(res)
}

## This version add one parameter maf_sids, which calcuates MAFs of SNP/Haps in a subset of sampels, e.g., BCa cases only
gdsCoxFit <- function(set.idx, gdsfn, pheno, yy, xx, type, frq.min, frq.max, maf_sids, haps, saveGeno=FALSE){
  
  gds <- seqOpen(gdsfn, allow.duplicate=TRUE)
  sid <- intersect(seqGetData(gds, 'sample.id'), pheno$f.eid)
  seqSetFilter(gds, sample.id=sid, variant.id=set.idx)
  pheno <- pheno[match(sid, pheno$f.eid), ]
  
  if(type=='snp'){
    geno <- seqGetData(gds, '$dosage_alt')
    
    if(missing(maf_sids)){
           eaf <- colMeans(geno, na.rm=T)/2
    }else{
           sid.idx <- which(pheno$f.eid%in%maf_sids)
           eaf <- colMeans(geno[sid.idx, ], na.rm=T)/2
    }
    maf <- eaf
    maf[which(eaf>0.5)] <- 1-eaf[which(eaf>0.5)]
    vidx <- which(maf>=frq.min&maf<frq.max)
	if(length(vidx)>0){
        if(length(vidx)==1){
               fit <- coxFit(geno[, vidx], nz=1, pheno[, yy], pheno[, xx])
               fit <- t(data.frame(fit))
        }else{
               geno <- geno[, vidx]
               fit <- t(apply(geno, 2, coxFit, nz=1, pheno[, yy], pheno[, xx]))
        }
        res <- data.frame(vid=seqGetData(gds, 'variant.id')[vidx], ea=seqGetData(gds, '$alt')[vidx], 
                           nea=seqGetData(gds, '$ref')[vidx], eaf=eaf[vidx], fit)
        if(saveGeno) res <- list(res=res, geno=geno)
	}else{
	      print('No SNP was tested, please check your data!')
         res <- NULL
	}
   
  }
  
  if(type=='hap'){
    hh <- seqGetData(gds, 'genotype')
    h1 <- apply(hh[1,,], 1, paste, collapse="")
    h2 <- apply(hh[2,,], 1, paste, collapse="")
    h1 <- paste('h',h1,sep='')
    h2 <- paste('h',h2,sep='')
    hh <- c(h1, h2)
    if(missing(haps)) haps <- unique(hh)
    nn <- length(hh)
    
    if(missing(maf_sids)){
          eaf <- as.numeric(sapply(haps, function(x){sum(hh==x)/nn}))
    }else{
          sid.idx <- which(pheno$f.eid%in%maf_sids)
          hh2 <- hh[c(sid.idx, nrow(pheno)+sid.idx)]
          nn2 <- length(hh2)
          eaf <- as.numeric(sapply(haps, function(x){sum(hh2==x)/nn2}))
    }
    maf <- eaf
    maf[which(eaf>0.5)] <- 1-eaf[which(eaf>0.5)]
    hidx <- which(maf>=frq.min&maf<frq.max)
    if(length(hidx)>0){
       hinfo <- data.frame(hap=haps[hidx], eaf.case=eaf[hidx])
       hinfo$eaf <- as.numeric(sapply(as.vector(unlist(haps[hidx])), function(x){sum(hh==x)/nn}))
       hinfo <- hinfo[order(-hinfo$eaf), ]
       hinfo <- data.frame(hid=paste('h',1:nrow(hinfo), sep=''), hinfo)
       chr <- as.numeric(seqGetData(gds, 'chromosome'))
       pos <- as.numeric(seqGetData(gds, 'position'))

       if(nrow(hinfo)==1){
              hmat <- as.numeric(h1==hinfo$hap)+as.numeric(h2==hinfo$hap)
              fit  <- coxFit(hmat, nz=1, pheno[, yy], pheno[, xx])
              fit  <- t(data.frame(fit))
       }else{
               hmat <- do.call('cbind', lapply(hinfo$hap, function(x){as.numeric(h1==x)})) + do.call('cbind', lapply(hinfo$hap, function(x){as.numeric(h2==x)}))
	        colnames(hmat) <- hinfo$hid
               fit <- t(apply(hmat, 2, coxFit, nz=1, pheno[, yy], pheno[, xx])) 
       }
       res <- data.frame(block=paste(chr[1], set.idx[1], set.idx[length(set.idx)], sep='_'), 
                        chr=chr[1], start=min(pos), end=max(pos), tot=nrow(hinfo), hinfo, fit)
       if(saveGeno)  res <- list(res=res, geno=hmat)

    }else{
	       print('No Haplotype was tested, please check your data!')
           res <- NULL
    }
}
seqClose(gds)
return(res)
}