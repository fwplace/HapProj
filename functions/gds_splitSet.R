args <- commandArgs(trailingOnly = TRUE)

splitSet <- function(vid, w, s, type){
  set.list <- list()
  st <- 1; ed <- st+w-1; i <- 1
  
  while(ed<=length(vid)){
    set.list[[i]] <- vid[st:ed]
    if(ed==length(vid)) break
    i <- i+1
    st <- st+s
    ed <- st+w-1
  }
 if(ed>length(vid)){
     ed <- length(vid)
     if(type=='hap') st <- ed-w+1
     if(type=='snp') st <- max(set.list[[i-1]])+1
     set.list[[i]] <- vid[st:ed]
  }
  return(set.list)
}

fn     <- args[1]
w      <- as.numeric(args[2])
s      <- as.numeric(args[3])
type   <- args[4]
bsize  <- as.numeric(args[5])
jobSplit <- as.numeric(args[6])
outdir   <- args[7]

load(fn)
vinfo <- get(gsub('.rdata','',basename(fn)))
for(chrid in 1:22){
   tt <- vinfo[vinfo$chr==chrid, ]
   set.list <- splitSet(tt$vid, w=w, s=s, type=type)
   set.idx  <- split(1:length(set.list), ceiling(seq_along(set.list)/bsize))
   sets <- lapply(set.idx, function(x){set.list[x]})
   bb <- split(1:length(sets), ceiling(seq_along(sets)/jobSplit))
   for(i in 1:length(bb)){
       bsets <- sets[bb[[i]]]
       save(bsets, file=paste(outdir, '/chr', chrid,'.bset.', i, '.rdata', sep=''))
       rm(bsets)
       print(paste('chrid=', chrid, ', bsets=',i,sep=''))
   }
}
