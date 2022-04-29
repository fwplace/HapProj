
gdsGeno <- function(vid, gdsfn, sid, type, haps=NULL){
  
  gds <- seqOpen(gdsfn, allow.duplicate=TRUE)
  seqSetFilter(gds, sample.id=sid, variant.id=vid)
  sid <- seqGetData(gds, 'sample.id')
  
  if(type=='snp'){
    if(length(vid)==1){
        geno <- as.numeric(seqGetData(gds, '$dosage_alt'))
        names(geno) <- sid
    }else{
        geno <- seqGetData(gds, '$dosage_alt')
        rownames(geno) <- sid
    }

  }
  if(type=='hap'){
    hh <- seqGetData(gds, 'genotype')
    h1 <- apply(hh[1,,], 1, paste, collapse="")
    h2 <- apply(hh[2,,], 1, paste, collapse="")
    h1 <- paste('h',h1,sep='')
    h2 <- paste('h',h2,sep='')
    hh <- c(h1, h2)
    if(length(haps)==1){
               geno <- as.numeric(h1==haps)+as.numeric(h2==haps)
               names(geno) <- sid
       }else{
               geno <- do.call('cbind', lapply(haps, function(x){as.numeric(h1==x)}))+
                       do.call('cbind', lapply(haps, function(x){as.numeric(h2==x)}))
              rownames(geno) <- sid
       }
    }
seqClose(gds)
return(geno)
}
