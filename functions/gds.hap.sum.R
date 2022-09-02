args <- commandArgs(trailingOnly = TRUE)
resPath <- as.character(args[1])
outfn1   <- as.character(args[2])
outfn2   <- as.character(args[3])
fns <- list.files(path=resPath, pattern='*.res.rdata', full.names=T)
tt <- basename(fns)
vv <- data.frame(chr=as.numeric(sapply(tt, function(x){substr(strsplit(x, '[.]')[[1]][1],4,5)})),
                 set=as.numeric(sapply(tt, function(x){strsplit(x, '[.]')[[1]][3]})),fn=fns)
vv <- vv[order(vv$chr, vv$set), ]
for(i in 1:length(fns)){
    load(vv$fn[i])
    res <- res[!is.na(res$p), ]
    res$hid <- paste(res$block, res$hid, sep='.')
    res2 <- res[order(res$p), ]
    res2 <- res2[!duplicated(res2$block), ]
    res2 <- res2[order(res2$start), ]
    if(i==1){
        write.table(res,  file=outfn1, sep='\t', col.names=T, row.names=F, quote=F, append=TRUE)
        write.table(res2, file=outfn2, sep='\t', col.names=T, row.names=F, quote=F, append=TRUE)
    }else{
        write.table(res,  file=outfn1, sep='\t', col.names=F, row.names=F, quote=F, append=TRUE)
        write.table(res2, file=outfn2, sep='\t', col.names=F, row.names=F, quote=F, append=TRUE)
    }
    rm(res, res2)
    print(i)
}
