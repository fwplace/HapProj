args <- commandArgs(trailingOnly = TRUE)
library(parallel)

resPath    <- as.character(args[1])
resPattern <- as.character(args[2])
pcut       <- as.numeric(args[3])
out1       <- as.character(args[4])
out2       <- as.character(args[5])

fns <- list.files(path=resPath, pattern=resPattern, full.names=T)
df1 <- df2 <- list()
for(i in 1:length(fns)){
    load(fns[i])
    res <- res[!is.na(res$p), ]
    res$hid <- paste(res$block, res$hid, sep='.')
    if(sum(res$p<pcut)>0){
       res.top <- res[res$p<pcut, ]
    }else{
        res.top <- NULL
    }
    res <- res[order(res$p), ]
    res <- res[!duplicated(res$block), ]
    res <- res[order(res$start), ]
    df1[[i]] <- res
    df2[[i]] <- res.top
    rm(res, res.top)
    print(i)
}

## get lead haplotype for each haplotype window
df1       <- do.call('rbind', df1)
df1$chr   <- as.numeric(df1$chr)
df1$start <- as.numeric(df1$start)
df1$end   <-  as.numeric(df1$end)
df1       <- df1[order(df1$chr, df1$start), ]
assign(basename(out1), df1)
save(list=basename(out1), file=paste(out1, '.rdata', sep=''))

## get all haplotypes pasing pcut 
df2       <- do.call('rbind', df2)
df2$chr   <- as.numeric(df2$chr)
df2$start <- as.numeric(df2$start)
df2$end   <-  as.numeric(df2$end)
df2       <- df2[order(df2$chr, df2$start), ]
assign(basename(out2), df2)
save(list=basename(out2), file=paste(out2, '.rdata', sep=''))

