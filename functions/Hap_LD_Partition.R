# Haplotype partition by LD
setwd("~/fwork/ProjHap/breast/analysis/DRIVE2/rare13")
system('mkdir partition')

library(SeqArray)
library(gpart)

## load haplotype
load('rare13.obj.rdata')
gdsfile <- '~/fwork/ProjHap/data/ukb.bca/gds/ukb.bca.hap.xxx.gds'
phenfile <- '~/fwork/ProjHap/data/ukb.bca/gds/ukb.bca.pheno.txt'
pheno <- read.table(file=phenfile, header=T, sep='\t')

for(i in 1:13){
dat <- rare13.obj[[i]]$snpset
gdsfn <- gsub('xxx',paste0('chr', dat$chr[1]), gdsfile)
gds <- seqOpen(gdsfn, allow.duplicate=TRUE)
sid <- intersect(seqGetData(gds, 'sample.id'), pheno$sid)
seqSetFilter(gds, sample.id=sid, variant.id=dat$vid)
geno <- seqGetData(gds, 'genotype')
hap <- rbind(geno[1,,], geno[2,,])

## Partition via BigLD
gdose <- geno[1,,]+ geno[2,,]
colnames(gdose) <- dat$idx
snp <- data.frame(chrN=dat$chr[1], rsID=dat$idx, bp=dat$pos)
fit <- BigLD(geno=gdose, SNPinfo=snp, clstgap=1e6, LD='r2', MAFcut=0.01)
tt <- as.numeric(unlist(lapply(1:nrow(fit), function(x){xx=as.numeric(fit[x,]);xx[2]:xx[3]})))
tmp <- setdiff(1:nrow(dat), tt)
df <- rbind(fit[,2:3], data.frame(start.index=tmp, end.index=tmp))
df <- df[order(df$start.index),]
df$size <- df$end.index-df$start.index+1
df <- data.frame(bid=rare13.obj[[i]]$sum$bid[1], block=1:nrow(df), df)
rownames(df) <- df$block
## calculate Frq
frq <- NULL
for(j in 1:nrow(df)){
	vidx <- df$start.index[j]:df$end.index[j]
	if(length(vidx)>1){
	      haps <- apply(hap[,vidx],1, paste, collapse="")
	      frq[j]  <- sum(haps==paste(dat[vidx,]$allele, collapse=''))/length(haps)
	}else{
	      frq[j]  <- sum(hap[,vidx]==dat$allele[vidx])/nrow(hap)
	}
	print(j)
}
df$frq <- frq
write.table(df, file=paste0('partition/ukbb_bigLD_locus',rare13.obj[[i]]$sum$bid[1],'.txt'),
            sep='\t', col.names=T, row.names=F, quote=F)
rm(fit, df)
seqClose(gds)
}

## estimate LD blocks by PLINK
for(i in 1:13){
dat <- rare13.obj[[i]]$snpset
gdsfn <- gsub('xxx',paste0('chr', dat$chr[1]), gdsfile)
gds <- seqOpen(gdsfn, allow.duplicate=TRUE)
sid <- intersect(seqGetData(gds, 'sample.id'), pheno$sid)
seqSetFilter(gds, sample.id=sid, variant.id=dat$vid)
geno <- seqGetData(gds, 'genotype')
hap <- rbind(geno[1,,], geno[2,,])
ped <- matrix(NA, nrow=nrow(geno[1,,]),ncol=2*nrow(dat))
ped[,seq(1,2*nrow(dat),2)] <- geno[1,,]
ped[,seq(2,2*nrow(dat),2)] <- geno[2,,]
ped <- ped + 1
ped <- data.frame(fid=0, iid=pheno$sid, pid=0, mid=0, sex=2, phen=pheno$diag+1, ped)
map <- data.frame(chr=dat$chr, snp=dat$idx, cm=0,  pos=dat$pos)
write.table(ped, file=paste0('partition/ukbb_locus',rare13.obj[[i]]$sum$bid[1],'.ped'),
            sep='\t', col.names=F, row.names=F, quote=F)
options(scipen=999)
write.table(map,file=paste0('partition/ukbb_locus',rare13.obj[[i]]$sum$bid[1],'.map'),
            sep='\t', col.names=F, row.names=F, quote=F)

system(paste0('plink --file ',paste0('partition/ukbb_locus',rare13.obj[[i]]$sum$bid[1]),
              ' --blocks --blocks-max-kb 1000 --blocks-min-maf 0.01 --out ',
              paste0('partition/ukbb_locus',rare13.obj[[i]]$sum$bid[1])))
tt <- read.table(file=paste0('partition/ukbb_locus', rare13.obj[[i]]$sum$bid[1],'.blocks.det'),header=T)
tt <- lapply(as.vector(unlist(tt$SNPS)), function(x){as.numeric(strsplit(x,'[|]')[[1]])})
df <- data.frame(start.index=sapply(tt, min), end.index=sapply(tt,max))
tt <- as.numeric(unlist(lapply(1:nrow(df), function(x){xx=as.numeric(df[x,]);xx[1]:xx[2]})))
tmp <- setdiff(1:nrow(dat), as.numeric(unlist(tt)))
df <- rbind(df, data.frame(start.index=tmp, end.index=tmp))
df <- df[order(df$start.index),]
df$size <- df$end.index-df$start.index+1
df <- data.frame(bid=rare13.obj[[i]]$sum$bid[1], block=1:nrow(df), df)
rownames(df) <- df$block
frq <- NULL
for(j in 1:nrow(df)){
  vidx <- df$start.index[j]:df$end.index[j]
  if(length(vidx)>1){
    haps <- apply(hap[,vidx],1, paste, collapse="")
    frq[j]  <- sum(haps==paste(dat[vidx,]$allele, collapse=''))/length(haps)
  }else{
    frq[j]  <- sum(hap[,vidx]==dat$allele[vidx])/nrow(hap)
  }
  print(j)
}
df$frq <- frq
write.table(df, file=paste0('partition/ukbb_PLINK_locus',rare13.obj[[i]]$sum$bid[1],'.txt'),
            sep='\t', col.names=T, row.names=F, quote=F)
seqClose(gds)
}

cd /Volumes/clusterhome/fwang2/ProjHap/breast/DRIVE2/rare13
R 

load('rare13.obj.rdata')

df <- list()
for(i in 1:13){
    df[[i]] <- read.table(file=paste0('partition/ukbb_bigLD_locus',
                                      rare13.obj[[i]]$sum$bid[1],'.txt'), header=T)
    print(i)
}
df <- do.call('rbind', df)
rare13_ukbb_bigLD <- df
save(rare13_ukbb_bigLD, file='partition/rare13_ukbb_bigLD.rdata')

df <- list()
for(i in 1:13){
    df[[i]] <- read.table(file=paste0('partition/ukbb_PLINK_locus',
                                      rare13.obj[[i]]$sum$bid[1],'.txt'), header=T)
    print(i)
}
df <- do.call('rbind', df)
rare13_ukbb_PLINK <- df
save(rare13_ukbb_PLINK, file='partition/rare13_ukbb_PLINK.rdata')

## plottting
setwd("~/fwork/ProjHap/breast/analysis/DRIVE2/rare13")
library(reshape)
library(ggplot2)
library(cowplot)
load('rare13.obj.rdata')
load('partition/rare13_ukbb_bigLD.rdata')
load('partition/rare13_ukbb_PLINK.rdata')

ld_plots <- list()
for(idx in 1:13){
  dat <- rare13.obj[[idx]]$snpset
  m1 <- as.matrix(rare13.obj[[idx]]$ukb.ld$LD.all)
  colnames(m1) <- rownames(m1) <- 1:nrow(dat)
  m2 <- melt(as.matrix(m1))
  res <- rare13.obj[[idx]]$sum
    ld_plots[[idx]] <- ggplot(m2, aes(X1, X2)) + geom_tile(aes(fill = value)) +
      scale_fill_gradient(low = "white", high = "red")+
      ggtitle(paste0('Locus ',res$bid[1],': ', res$block[1],' (W=',res$run[1],')')) +
      guides(fill=guide_legend(title=expression(D^"'"/~italic(r)^2)))+
      geom_segment(x=1, y = 1, xend = nrow(dat), yend = nrow(dat), color='black', alpha=0.5, size=0.6)+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.border = element_blank(), panel.background = element_blank(),
                         axis.line = element_blank(),plot.title = element_text(size = 8, face='bold',margin=margin(0,0,0,0)),
                         axis.text = element_blank(), axis.ticks = element_blank(), 
                         axis.title = element_blank(),
                         plot.margin = unit(c(0.1,0.1,0,0.1), "cm"))
}

ld_bigLD <- list()
for(idx in 1:13){
res <- rare13.obj[[idx]]$sum
dat <- rare13.obj[[idx]]$snp
df <- rare13_ukbb_bigLD[rare13_ukbb_bigLD$bid==res$bid[1], ]
pp <- ggplot(data=df[df$size==1, ], aes(x=start.index, y=log10(frq))) + 
  geom_point(size=0.5, colour='blue')+
    scale_x_continuous(name='SNP Index', limits = c(0, nrow(dat)))+
    scale_y_continuous(name='Frq', limits = c(-4, 0), 
                       breaks = log10(c(10^-4, 10^-3, 10^-2, 0.1, 1)),
                       labels = c(expression(10^-4), expression(10^-3), expression(10^-2), 0.1, 1.0))+
  geom_segment(data=df[df$size>1, ], 
               aes(x=start.index , y = log10(frq), xend = end.index, yend = log10(frq)), 
               color='darkcyan', alpha=0.8, size=1.5) +
  geom_vline(xintercept=dat[dat$validation==1, ]$idx, color = "gray",size=0.4)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), panel.background = element_blank(),
                     axis.line = element_line(colour = "black", size=0.5),
                     axis.text = element_text(size=8,colour = "black"), axis.title = element_text(size=8),
                     axis.ticks.length = unit(0.05,"cm"), plot.margin = unit(c(0,0.1,0.1,0.1), "cm"))
 pp1 <- ld_plots[[idx]] + theme(legend.position='none')
 ld_bigLD[[idx]] <- plot_grid(pp1, pp, ncol=1, align = 'v')
}

png(file='rare13_LDplots_bigLD.png',height=8,width=7, unit='in', res=600)
legend_pp <- get_legend(ld_plots[[1]])
plot_grid(ld_bigLD[[1]], ld_bigLD[[2]],ld_bigLD[[3]],
          ld_bigLD[[4]], ld_bigLD[[5]],ld_bigLD[[6]],
          ld_bigLD[[7]], ld_bigLD[[8]],ld_bigLD[[9]],
          ld_bigLD[[10]], ld_bigLD[[11]],ld_bigLD[[12]],
          ld_bigLD[[13]], legend_pp, ncol=3, labels=c(LETTERS[1:13],''), label_size = 10)
dev.off()

pdf(file='rare13_LDplots_bigLD.pdf',height=8,width=7)
plot_grid(ld_bigLD[[1]], ld_bigLD[[2]],ld_bigLD[[3]],
          ld_bigLD[[4]], ld_bigLD[[5]],ld_bigLD[[6]],
          ld_bigLD[[7]], ld_bigLD[[8]],ld_bigLD[[9]],
          ld_bigLD[[10]], ld_bigLD[[11]],ld_bigLD[[12]],
          ld_bigLD[[13]], ncol=3, labels=LETTERS[1:13], label_size = 10)
dev.off()


ld_PLINK <- list()
for(idx in 1:13){
  res <- rare13.obj[[idx]]$sum
  dat <- rare13.obj[[idx]]$snpset
  df <- rare13_ukbb_PLINK[rare13_ukbb_PLINK$bid==res$bid[1], ]
pp <- ggplot(data=df[df$size==1, ], aes(x=start.index, y=log10(frq))) + 
  geom_point(size=0.5, colour='blue')+
    scale_x_continuous(name='SNP Index', limits = c(0, nrow(dat)))+
    scale_y_continuous(name='Frq', limits = c(-4, 0), 
                       breaks = log10(c(10^-4, 10^-3, 10^-2, 0.1, 1)),
                       labels = c(expression(10^-4), expression(10^-3), expression(10^-2), 0.1, 1.0))+
  geom_segment(data=df[df$size>1, ], 
               aes(x=start.index , y = log10(frq), xend = end.index, yend = log10(frq)), 
               color='darkcyan', alpha=0.8, size=1.5) +
  geom_vline(xintercept=dat[dat$validation==1, ]$idx, color = "gray",size=0.4)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), panel.background = element_blank(),
                     axis.line = element_line(colour = "black", size=0.5),
                     axis.text = element_text(size=8,colour = "black"), axis.title = element_text(size=8),
                     axis.ticks.length = unit(0.05,"cm"), plot.margin = unit(c(0,0.1,0.1,0.1), "cm"))
 pp1 <- ld_plots[[idx]] + theme(legend.position='none')
 ld_PLINK[[idx]] <- plot_grid(pp1, pp, ncol=1, align = 'v')
}


png(file='rare13_LDplots_PLINK.png',height=8,width=7, unit='in', res=600)
legend_pp <- get_legend(ld_plots[[1]])
plot_grid(ld_PLINK[[1]], ld_PLINK[[2]],ld_PLINK[[3]],
          ld_PLINK[[4]], ld_PLINK[[5]],ld_PLINK[[6]],
          ld_PLINK[[7]], ld_PLINK[[8]],ld_PLINK[[9]],
          ld_PLINK[[10]], ld_PLINK[[11]],ld_PLINK[[12]],
          ld_PLINK[[13]], legend_pp, ncol=3, labels=LETTERS[1:13], label_size = 10)
dev.off()

pdf(file='rare13_LDplots_PLINK.pdf',height=8,width=7)
plot_grid(ld_PLINK[[1]], ld_PLINK[[2]],ld_PLINK[[3]],
          ld_PLINK[[4]], ld_PLINK[[5]],ld_PLINK[[6]],
          ld_PLINK[[7]], ld_PLINK[[8]],ld_PLINK[[9]],
          ld_PLINK[[10]], ld_PLINK[[11]],ld_PLINK[[12]],
          ld_PLINK[[13]], ncol=3, labels=LETTERS[1:13], label_size = 10)
dev.off()




