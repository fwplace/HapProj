ld_find <- function(dat, fn.ld){
    tt <- read.table(file=fn.ld, header=T)[, c(1:5,12)]
    tt <- tt[order(tt$CHR, tt$BP), ]
    sp2 <- list()
    for(j in 1:nrow(tt)){
          if(tt$SP2[j]=='NONE') sp2[[j]] <- tt$SNP[j]
          if(tt$SP2[j]!='NONE'){
               tmp <- strsplit(tt$SP2[[j]], '[,]')[[1]]
               tmp <- as.vector(unlist(sapply(tmp, function(x){strsplit(x,'[()]')[[1]][1]})))
               sp2[[j]] <- c(tt$SNP[j], tmp)
          }
     }
     ld <- data.frame(ld=rep(1:nrow(tt), lapply(sp2,length)),
                      sp1=rep(tt$SNP,lapply(sp2,length)),sp2=as.vector(unlist(sp2)))
     ld <- ld[match(dat$hid, ld$sp2), ]
     dat$ld <- ld$ld
     tt <- dat[order(dat$pros.p), ]
     tt <- tt[!duplicated(tt$ld), ]
     dat$lead <- as.numeric(dat$hid%in%tt$hid)
    return(dat)
}

