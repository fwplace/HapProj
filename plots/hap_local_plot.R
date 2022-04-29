## generating local association plots

## 1. prepare data tracks
cd ~/fwork/ProjHap/breast/rare13
mkdir tracks

## SNP-GWAS
load('~/fwork/ProjHap/breast/analysis/ukb.bca.typed.rdata')
snp.gwas <- ukb.bca.typed[, c(1:6,28:31)]
colnames(snp.gwas)[7:10] <- c('eaf','hr','se','p')
save(snp.gwas, file='tracks/snp.gwas.rdata')
rm(ukb.bca.typed)

## genetic maps
load('~/fwork/PublicData/Bherer_genetic_map/Bherer_female_recomb.rdata')
recomb.female <- Bherer_female_recomb
save(recomb.female, file='tracks/recomb.female.rdata')
rm(Bherer_female_recomb)

## TAD
load('~/fwork/PublicData/tad/tads.hg19.rdata')
tad.HMEC <- tads.hg19[tads.hg19$cell=='HMEC_Lieberman', ]
save(tad.HMEC, file='tracks/tad.HMEC.rdata')
rm(tads.hg19)

## chromatin loops
load('~/fwork/PublicData/tad/loops.hg19.rdata')
tt <- loops.hg19[loops.hg19$cell%in%c('Rao_2014.HMEC','ENCODE3.T47D'), ]
tt$name <- paste(tt$chr,tt$st1,tt$ed1,tt$st2,tt$ed2,sep=':')
tt <- tt[!duplicated(tt$name), ]
tt <- tt[order(tt$chr,tt$st1), ]
loops.breast <- tt[, c(1,3:7)]
save(loops.breast, file='tracks/loops.breast.rdata')
rm(tt, loops.hg19)

## gene track within +/- 200Kb 
load('~/fwork/PublicData/gencode.v37.hg19.rdata')
load('anno2/vinfo.rdata')
hdat <- vinfo[!duplicated(vinfo$bid), c(1:14,21)]
df <- list()
for(i in 1:13){
    tt <- gencode.v37.hg19[gencode.v37.hg19$seqnames==hdat$chr[i], ]
    d1 <- tt$start-(hdat$start[i]-2e5)
    d2 <- tt$start-(hdat$end[i]+2e5)
    d3 <- tt$end-(hdat$start[i]-2e5)
    d4 <- tt$end-(hdat$end[i]+2e5)
    tmp <- tt[sign(d1)!=sign(d2)|sign(d3)!=sign(d4), ]
    if(nrow(tmp)>0){
        df[[i]] <- data.frame(bid=hdat$bid[i], tmp[order(tmp$start), ])
    }else{
        df[[i]] <- NULL
    }
print(i)
}
df <- do.call('rbind', df)
write.csv(df, file='tracks/rare13_genes_200k.csv')

## plot...
library(ggplot2)
library(cowplot)
library(ggrepel)
load('anno2/vinfo.rdata')
load('tracks/snp.gwas.rdata')
load('tracks/recomb.female.rdata')
load('tracks/tad.HMEC.rdata')
load('tracks/loops.breast.rdata')
rare13.genes <- read.table(file='tracks/rare13_genes_200k.txt',header=T)
hdat <- vinfo[!duplicated(vinfo$bid), c(1:21,32)]
hdat <- hdat[order(hdat$bid), ]

## bid
assoc_plots <- list()
dist0 <- 2e5
for(i in 1:13){
chr <- hdat$chr[i]
st  <- hdat$start[i]-dist0
ed  <- hdat$end[i]+dist0
d0 <- vinfo[vinfo$bid==hdat$bid[i], ]
d1 <- snp.gwas[snp.gwas$chr==chr&snp.gwas$pos>=st&snp.gwas$pos<=ed, ]
d2 <- recomb.female[recomb.female$chr==chr&recomb.female$pos>=st&recomb.female$pos<=ed, ]
d2$rate[d2$rate>100] <- 100
d3 <- tad.HMEC[tad.HMEC$chr==chr, ]
tmp1 <- d3$start<=st&d3$end>=st|d3$start>=st&d3$start<=ed
d3 <- d3[tmp1, ]
d3$start[d3$start<st] <- st
d3$end[d3$end>ed] <- ed
d4 <- loops.breast[loops.breast$chr==chr, ]
tmp1 <- d4$st1<=st&d4$ed1>=st|d4$st1>=st&d4$st1<=ed
tmp2 <- d4$st2<=st&d4$ed2>=st|d4$st2>=st&d4$st2<=ed
d4 <- d4[tmp1|tmp2, ]
d4$st1[d4$st1<st] <- st
d4$ed1[d4$ed1>ed] <- ed
d4$st2[d4$st2<st] <- st
d4$ed2[d4$ed2>ed] <- ed
d5 <- rare13.genes[rare13.genes$bid==hdat$bid[i], ]
d5$start[d5$start<st] <- st
d5$end[d5$end>ed] <- ed
ymax <- 16

pp <- ggplot(data=d1, aes(x=pos/1e6, y=-log10(p))) + 
       geom_point(size=0.6, colour='black') +
       scale_x_continuous(name='Mb', limits = c(st/1e6, ed/1e6))+
       scale_y_continuous(name='-Log10(P)', limits = c(-8, ymax), breaks=seq(0,ymax,4), 
                              sec.axis=sec_axis(~.*(100/ymax), breaks=seq(0,100,20), name='cM/Mb'))+
       geom_line(data=d2, aes(x=pos/1e6, y=rate*ymax/100), color="brown", size=0.4, alpha=0.6)+
       ggtitle(paste0('Bid=',hdat$bid[i],': ', hdat$block[i],'(',hdat$cytoband[i],')')) +
       geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "black", alpha=0.5,size=0.4)+
       geom_segment(data=data.frame(x1=min(d0$pos), x2=max(d0$pos), y=-log10(hdat$comb.cox[i])), 
                 aes(x=x1/1e6 , y = y, xend = x2/1e6, yend = y), color='orange', size=1.5) +
       geom_segment(data=d0, aes(x=pos/1e6, y = -log10(comb.cox)-0.3, xend = pos/1e6, yend = -log10(comb.cox)+0.3), 
                 color='black', alpha=0.5, size=0.3) +
       geom_segment(data=data.frame(x1=min(d0$pos[d0$validation==1]), x2=max(d0$pos[d0$validation==1]), y=-log10(hdat$ukb.cox[i])), 
                 aes(x=x1/1e6 , y = y, xend = x2/1e6, yend = y), color='purple', size=1.5)+
       geom_segment(data=d0[d0$validation==1, ], aes(x=pos/1e6, y = -log10(ukb.cox)-0.3, xend = pos/1e6, yend = -log10(ukb.cox)+0.3), 
                 color='black', alpha=0.5, size=0.5)

if(nrow(d3)>0) pp <- pp + geom_segment(data=d3, aes(x=start/1e6 , y = -1, xend = end/1e6, yend = -1), color='darkgreen', alpha=0.8, size=1)

if(nrow(d4)>0) pp <- pp + geom_segment(data=d4, aes(x=st1/1e6 , y = -2, xend = ed1/1e6, yend = -2), color='magenta', size=2)+
                          geom_segment(data=d4, aes(x=st2/1e6 , y = -2, xend = ed2/1e6, yend = -2), color='magenta', size=2)+
                          geom_curve(data=d4, aes(x=(st1+ed1)/2e6 , y = -2.1, xend = (st2+ed2)/2e6, yend = -2.1), color='magenta', curvature = 0.5, size=0.3)

if(sum(d5$strand=='+')>0) pp <- pp + geom_segment(data=d5[d5$strand=='+', ], aes(x=start/1e6, y=group-8, xend=end/1e6, yend=group-8), 
                          arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.8)

if(sum(d5$strand=='-')>0) pp <- pp + geom_segment(data=d5[d5$strand=='-', ], aes(x=end/1e6, y=group-8, xend=start/1e6, yend=group-8), 
                          arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.8)

if(nrow(d5)>0) pp <- pp + geom_text_repel(data=d5, aes(x=(start+end)/2e6, y=group-8, label=gene_name), size=2, 
                                          box.padding = unit(0.05, "lines"), point.padding = unit(0.05, "lines"))

assoc_plots[[i]] <- pp + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.border = element_blank(), panel.background = element_blank(),
                          axis.line = element_line(colour = "black"), plot.title = element_text(size = 10),
                          axis.text = element_text(size=10,colour = "black"), axis.title = element_text(size=10),
                          axis.title.x = element_text(vjust=0.6),axis.title.y = element_text(vjust=0.6),
                          axis.ticks.length = unit(0.25,"cm"), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
}
pdf(file='Rare13_assoc_plots.pdf',height=12,width=10)
plot_grid(assoc_plots[[1]], assoc_plots[[2]], assoc_plots[[3]],
          assoc_plots[[4]],assoc_plots[[5]],assoc_plots[[6]],
          assoc_plots[[7]],assoc_plots[[8]],assoc_plots[[9]],
          assoc_plots[[10]],assoc_plots[[11]],assoc_plots[[12]],
          assoc_plots[[13]], ncol=3, labels=LETTERS[1:13],align = 'hv')
dev.off()

save(assoc_plots, file='tracks/assoc_plots.rdata')



## single locus
dist0 <- 2e5
i=8
chr <- hdat$chr[i]
st  <- hdat$start[i]-dist0
ed  <- hdat$end[i]+dist0
d0 <- vinfo[vinfo$bid==hdat$bid[i], ]
d1 <- snp.gwas[snp.gwas$chr==chr&snp.gwas$pos>=st&snp.gwas$pos<=ed, ]
d2 <- recomb.female[recomb.female$chr==chr&recomb.female$pos>=st&recomb.female$pos<=ed, ]
d2$rate[d2$rate>100] <- 100
d3 <- tad.HMEC[tad.HMEC$chr==chr, ]
tmp1 <- d3$start<=st&d3$end>=st|d3$start>=st&d3$start<=ed
d3 <- d3[tmp1, ]
d3$start[d3$start<st] <- st
d3$end[d3$end>ed] <- ed
d4 <- loops.breast[loops.breast$chr==chr, ]
tmp1 <- d4$st1<=st&d4$ed1>=st|d4$st1>=st&d4$st1<=ed
tmp2 <- d4$st2<=st&d4$ed2>=st|d4$st2>=st&d4$st2<=ed
d4 <- d4[tmp1|tmp2, ]
d4$st1[d4$st1<st] <- st
d4$ed1[d4$ed1>ed] <- ed
d4$st2[d4$st2<st] <- st
d4$ed2[d4$ed2>ed] <- ed
d5 <- rare13.genes[rare13.genes$bid==hdat$bid[i], ]
d5$start[d5$start<st] <- st
d5$end[d5$end>ed] <- ed
ymax <- 16
pp <- ggplot(data=d1, aes(x=pos/1e6, y=-log10(p))) + 
        geom_point(size=1, colour='black') +
        scale_x_continuous(name='Mb', limits = c(st/1e6, ed/1e6))+
        scale_y_continuous(name='-Log10(P)', limits = c(-8, ymax), breaks=seq(0,ymax,4), 
                           sec.axis=sec_axis(~.*(100/ymax), breaks=seq(0,100,20), name='cM/Mb'))+
        geom_line(data=d2, aes(x=pos/1e6, y=rate*ymax/100), color="brown", size=1, alpha=0.6)+
        ggtitle(paste0('Bid=',hdat$bid[i],': ', hdat$block[i],'(',hdat$cytoband[i],')')) +
        geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "black", alpha=0.5,size=0.8)+
        geom_segment(data=data.frame(x1=min(d0$pos), x2=max(d0$pos), y=-log10(hdat$comb.cox[i])), 
                     aes(x=x1/1e6 , y = y, xend = x2/1e6, yend = y), color='orange', size=1) +
        geom_segment(data=d0, aes(x=pos/1e6, y = -log10(comb.cox)-0.3, xend = pos/1e6, yend = -log10(comb.cox)+0.3), 
                     color='black', alpha=0.5, size=1) +
        geom_segment(data=data.frame(x1=min(d0$pos[d0$validation==1]), x2=max(d0$pos[d0$validation==1]), y=-log10(hdat$ukb.cox[i])), 
                     aes(x=x1/1e6 , y = y, xend = x2/1e6, yend = y), color='purple', size=1)+
        geom_segment(data=d0[d0$validation==1, ], aes(x=pos/1e6, y = -log10(ukb.cox)-0.3, xend = pos/1e6, yend = -log10(ukb.cox)+0.3), 
                     color='black', alpha=0.5, size=1)
    
if(nrow(d3)>0) pp <- pp + geom_segment(data=d3, aes(x=start/1e6 , y = -1, xend = end/1e6, yend = -1), color='darkgreen', alpha=0.8, size=1)
    
if(nrow(d4)>0) pp <- pp + geom_segment(data=d4, aes(x=st1/1e6 , y = -2, xend = ed1/1e6, yend = -2), color='magenta', size=2)+
        geom_segment(data=d4, aes(x=st2/1e6 , y = -2, xend = ed2/1e6, yend = -2), color='magenta', size=2)+
        geom_curve(data=d4, aes(x=(st1+ed1)/2e6 , y = -2.1, xend = (st2+ed2)/2e6, yend = -2.1), color='magenta', curvature = 0.5, size=0.3)
    
if(sum(d5$strand=='+')>0) pp <- pp + geom_segment(data=d5[d5$strand=='+', ], aes(x=start/1e6, y=group-8, xend=end/1e6, yend=group-8), 
                                                      arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.8)
    
if(sum(d5$strand=='-')>0) pp <- pp + geom_segment(data=d5[d5$strand=='-', ], aes(x=end/1e6, y=group-8, xend=start/1e6, yend=group-8), 
                                                      arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.8)
    
if(nrow(d5)>0) pp <- pp + geom_text_repel(data=d5, aes(x=(start+end)/2e6, y=group-8, label=gene_name), size=3, 
                                              box.padding = unit(0.05, "lines"), point.padding = unit(0.05, "lines"))

pp <- pp + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                panel.border = element_blank(), panel.background = element_blank(),
                                                axis.line = element_line(colour = "black"), plot.title = element_text(size = 14),
                                                axis.text = element_text(size=14,colour = "black"), axis.title = element_text(size=14),
                                                axis.title.x = element_text(vjust=0.6),axis.title.y = element_text(vjust=0.6),
                                                axis.ticks.length = unit(0.25,"cm"), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

pdf(file='Rare13_assoc_plots_bid15.pdf',height=4,width=6)
pp
dev.off()


## summary
load('tracks/assoc_plots.rdata')
anno <- read.csv(file='rare13_annotation_summary.csv')[, -1]

df <- data.frame(x=rep(1:13,5), name=rep(anno$bid, 5),
                 type=rep(c('eQTL','Enhancer','TADBoundary','Chromatin Loop','Motif Change'), each=13), 
                 frq=c(anno$hits.gtex8.bm/anno$w, anno$hits.enhancer.breast/anno$w, anno$hits.tad.boundary.20k/anno$w, 
                   anno$hits.loop/anno$w, anno$hits.motif.diff3/anno$w))
df$type <- factor(df$type, levels=c('eQTL','Enhancer','TADBoundary','Chromatin Loop','Motif Change'))


pp <- ggplot(data=df, aes(x=x, y=frq, col=type)) + 
    geom_line(size=0.6) +
    scale_x_continuous(name='Haplotype Locus (Bid)', breaks=seq(1,13,1), labels = anno$bid)+
    scale_y_continuous(name='% Annotated SNPs', limits=c(0,1), breaks=seq(0,1,0.2))+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.border = element_blank(), panel.background = element_blank(),
                       axis.line = element_line(colour = "black"), plot.title = element_text(size = 10),
                       axis.text = element_text(size=10,colour = "black"), axis.title = element_text(size=10),
                       axis.title.x = element_text(vjust=0.6),axis.title.y = element_text(vjust=0.6),
                       axis.ticks.length = unit(0.25,"cm"), legend.title = element_blank()) 
    
pdf(file='Rare13_assoc_plots_v2.pdf',height=12,width=10)

top <- plot_grid(assoc_plots[[1]], assoc_plots[[2]], assoc_plots[[3]],
          assoc_plots[[4]],assoc_plots[[5]],assoc_plots[[6]],
          assoc_plots[[7]],assoc_plots[[8]],assoc_plots[[9]],
          assoc_plots[[10]],assoc_plots[[11]],assoc_plots[[12]], ncol=3, 
          labels=LETTERS[1:12],align = 'hv')
bottom <- plot_grid(assoc_plots[[13]], pp, rel_widths = c(1,2), nrow=1, labels=LETTERS[13:14],align = 'hv')
plot_grid(top, bottom,rel_heights = c(4,1),ncol=1)
dev.off()

