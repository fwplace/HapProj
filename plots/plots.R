setwd('~/fwork/ProjHap/breast')
library(ggplot2)

## Figure S1
cd ~/fwork/ProjHap/breast
R

library(ggplot2)
library(cowplot)

load('~/fwork/ProjHap/data/ukb.bca/ukb.bca.1kg.sinfo.rdata')
ukb.bca.1kg.sinfo$super_pop <- gsub('Study', 'UKBB', ukb.bca.1kg.sinfo$super_pop)
ukb.bca.1kg.sinfo$super_pop <- factor(ukb.bca.1kg.sinfo$super_pop, levels=c('AFR','AMR','EAS','EUR','SAS','UKBB'))
p1 <- ggplot(ukb.bca.1kg.sinfo, aes(x=pc1, y=pc2, color=as.factor(super_pop))) +
      geom_point(size=0.7)+
      scale_colour_discrete('')+
      labs(x ="PC1", y = "PC2")+
      theme_classic()+
      theme(plot.title = element_text(size = 12),
            axis.text = element_text(size=12,colour = "black"),
            axis.title = element_text(size=12),
            axis.ticks.length = unit(0.25,"cm"),
            legend.title=element_text(size=10,colour = "black"), 
            legend.position=c(0.85,0.35),legend.spacing.x = unit(0.05, 'cm'),
            legend.text=element_text(size=10))

p2 <- ggplot(ukb.bca.1kg.sinfo, aes(x=pc1, y=pc3, color=as.factor(super_pop))) +
      geom_point(size=0.7)+
      scale_colour_discrete('')+
      labs(x ="PC1", y = "PC3")+
      theme_classic()+
      theme(plot.title = element_text(size = 12),
            axis.text = element_text(size=12,colour = "black"),
            axis.title = element_text(size=12),
            axis.ticks.length = unit(0.25,"cm"),
            legend.title=element_text(size=10,colour = "black"), 
            legend.position=c(0.85,0.35),legend.spacing.x = unit(0.05, 'cm'),
            legend.text=element_text(size=10))

pdf(file='summary/plots/UKBB_PCA.pdf', height=3.5,width=7)
plot_grid(p1, p2,  labels = c('A', 'B'))
dev.off()


## Figure S2
cd ~/fwork/ProjHap/breast
R

library(ggplot2)
library(cowplot)

load('~/fwork/ProjHap/hap/phasing/dbgap28544/plink/dbgap28544.1kg.sinfo.rdata')
dbgap28544.1kg.sinfo$pop <- gsub('Study', 'DRIVE', dbgap28544.1kg.sinfo$pop)
dbgap28544.1kg.sinfo$pop <- factor(dbgap28544.1kg.sinfo$pop, levels=c('AFR','AMR','EAS','EUR','SAS','DRIVE'))
rls   <- read.table(file='~/fwork/ProjHap/hap/phasing/dbgap28544/plink/dbgap28544.f59958.kin0.rls2', header=T)

p1 <- ggplot(dbgap28544.1kg.sinfo, aes(x=pc1, y=pc2, color=pop)) +
      geom_point(size=0.7)+
      scale_colour_discrete('')+
      labs(title='', x ="PC1", y = "PC2")+
      geom_vline(xintercept=0.0025, linetype="dashed", color="gray") +
      geom_hline(yintercept=-0.0070, linetype="dashed", color="gray") +
      theme_classic()+
      theme(plot.title = element_text(size = 12),
            axis.text = element_text(size=12,colour = "black"),
            axis.title = element_text(size=12),
            axis.ticks.length = unit(0.25,"cm"),
            legend.title=element_text(size=10,colour = "black"), 
            legend.position=c(0.85,0.35),legend.spacing.x = unit(0.03, 'cm'),
            legend.text=element_text(size=10))
            
p2 <- ggplot(rls, aes(IBS0, Kinship)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="gray") +
    geom_point(color='black',alpha=0.5)+
    labs(title='', ylab("kinship estimate"), xlab='IBS0')+
    theme_classic()+
          theme(plot.title = element_text(size = 12),
            axis.text = element_text(size=12,colour = "black"),
            axis.title = element_text(size=12),
            axis.ticks.length = unit(0.25,"cm"))

pdf(file='summary/plots/Drive_validation_qc.pdf',height=3.5,width=7)
plot_grid(p2, p1, labels=c('A','B'))
dev.off()

### Figure S3
library(ggplot2)
library(cowplot)
load('ukb.bca.rare13.rdata')
df <- do.call('rbind', ukb.bca.rare13$hap2snp)
df <- df[!duplicated(df$snp), ] ## 1813 SNPs
df <- df[!is.na(df$imp.rsq), ]
df$comb.maf <- df$comb.af
df$comb.maf[which(df$comb.af>0.50)] <- 1-df$comb.af[which(df$comb.af>0.50)]
df$grp <- ifelse(df$imp.maf>=0.01, 'common(n=1,390)', 'rare(n=126)')
df$grp <- factor(df$grp, levels=c('common(n=1,390)', 'rare(n=126)'))

p1 <- ggplot(df, aes(x=imp.maf, y=imp.rsq, color=grp)) + geom_point(size=0.7) +
    labs(title='DRIVE TOPMED-Imputation', x='MAF', y='r2')+
    scale_y_continuous(limits = c(0.2, 1.0))+
    geom_hline(yintercept=0.80, linetype="dashed", color="black") +
        theme_classic()+
          theme(plot.title = element_text(size = 10),
            axis.text = element_text(size=10,colour = "black"),
            axis.title = element_text(size=10),
            legend.position=c(0.8,0.2), legend.title=element_blank(),
            legend.key=element_blank(), legend.spacing.x = unit(0.0, 'cm'),
            axis.ticks.length = unit(0.25,"cm"))

p2 <- ggplot(df, aes(x=comb.maf, y=imp.maf)) + geom_point(size=0.7) +
    labs(title='1,516 SNPs on 193 Rare Haplotypes', x='DRIVE, MAF', y='UKBB, MAF')+
    geom_abline(intercept = 0, slope = 1, color='red',size=0.5) +
        theme_classic()+
          theme(plot.title = element_text(size = 10),
            axis.text = element_text(size=10,colour = "black"),
            axis.title = element_text(size=10),legend.position='none',
            axis.ticks.length = unit(0.25,"cm"))

pdf(file='summary/plots/SNP_MAFvsr2.pdf',height=3.5,width=7)
plot_grid(p1, p2, labels=c('A','B'))
dev.off()

## Figure 1B
load('validation/res/ukb.bca.top436.rdata')
res <-  ukb.bca.top436
w <- c(5,10,20,30,50,100,250,500)
t1 <- as.numeric(sapply(w, function(x){sum(res$run==x&res$retro.afcase>100)}))
t2 <- as.numeric(sapply(w, function(x){sum(res$run==x&res$retro.afcase<=100)}))
df <- data.frame(run=rep(w,2),type=rep(c('common','rare'),each=8), count=c(t1,t2))
df$run <- factor(df$run, levels=w)
df$type <- factor(df$type, levels=c('common','rare'))
df2 <- data.frame(run=w, tot=df$count[1:8]+df$count[9:16])
df2$run <- factor(df2$run, levels=w)
ggplot(df, aes(fill=type, y=count, x=run, label=count)) + 
    geom_bar(position="stack", stat="identity")+ 
    geom_text(data=df2, aes(x=run, y=tot, label=tot, fill=NULL), vjust = -0.3, size = 4)+
    xlab('W') + scale_y_continuous(name='# Of Risk Haplotypes ', limits = c(0, 150))+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), panel.border = element_blank(),
                        axis.line = element_line(color='black'),
                        plot.title = element_text(size = 12),
                        axis.text = element_text(size=12,colour = "black"),
                        axis.title = element_text(size=12),
                        axis.ticks.length = unit(0.25,"cm"),
                        legend.title=element_blank(), 
                        legend.position=c(0.4,0.9), 
                        legend.spacing.x = unit(0.05, 'cm'),
                        legend.text=element_text(size=10))


## LD heatmap plots in UKBB

## hap-assoc plots
library(ggplot2)
library(cowplot)
library(ggrepel)
load('rare11/rare11.obj2.rdata')
load('ukb.bca.typed.rdata')
load('anno/gencode.v37.hg19.rdata')
load('anno/Bherer_female_recomb.rdata')
rare11.gene.200k <- read.table(file='rare11/rare11_genes_200k.txt',sep='\t',header=T)

pp <- list()
for(idx in 1:11){
d1 <- rare11.obj2[[idx]]$sum
d2 <- rare11.obj2[[idx]]$snp.anno
d3 <- rare11.obj2[[idx]]$val
dist0 <- 2e5
chr <- d1$chr
st  <- d1$start-dist0
ed  <- d1$end+dist0
d4  <- Bherer_female_recomb[Bherer_female_recomb$chr==chr&Bherer_female_recomb$pos>=st&Bherer_female_recomb$pos<=ed, ]
d4$rate[d4$rate>100] <- 100
d5  <- ukb.bca.typed[ukb.bca.typed$chr==chr&ukb.bca.typed$pos>=st&ukb.bca.typed$pos<=ed, ]
ymax <- round(-log10(d1$comb.cox),digits=0)+1
d6 <- rare11.gene.200k[rare11.gene.200k$bid==rare11.obj2[[idx]]$sum$bid, ]
d6$start[d6$start<st] <- st
d6$end[d6$end>ed] <- ed

# png(file=paste0('summary/plots/hap_assoc_bid',d1$bid[1],'.png'),height=3,width=6,unit='in', res=300)
p1 <- ggplot(data=d5, aes(x=pos/1e6, y=-log10(comb.cox))) +
      geom_point(size=0.8, colour='black', alpha=0.8) +
      scale_x_continuous(name=paste0('Chr ',chr, ' (Mb)'), limits = c(st/1e6, ed/1e6))+
      scale_y_continuous(name='-Log10(P)', limits = c(-4, ymax), breaks=seq(0,ymax,5), 
                         sec.axis=sec_axis(~.*(100/ymax), breaks=seq(0,100,20), name='Recomb Rate, cM/Mb'))+
      geom_line(data=d4, aes(x=pos/1e6, y=rate*ymax/100), color="brown", size=0.4)+
      ggtitle(paste0('Bid=',d1$bid,': ',d1$block)) +
      geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "black", alpha=0.5,size=0.4)+
      geom_segment(aes(x =min(d2$pos)/1e6 , y = -log10(d1$comb.cox), xend = max(d2$pos)/1e6, yend = -log10(d1$comb.cox)), color='orange', size=1) +
      geom_segment(data=d2, aes(x=pos/1e6, y = -log10(d1$comb.cox)-0.3, xend = pos/1e6, yend = -log10(d1$comb.cox)+0.3), color='black', alpha=0.8, size=0.3) +
      geom_segment(aes(x =min(d2$pos[d2$validation==1])/1e6 , y = -log10(d1$ukb.cox), 
                   xend = max(d2$pos[d2$validation==1])/1e6, yend = -log10(d1$ukb.cox)), color='purple', size=1) +
      geom_segment(data=d2[d2$validation==1, ], aes(x=pos/1e6, y = -log10(d1$ukb.cox)-0.3, 
                  xend = pos/1e6, yend = -log10(d1$ukb.cox)+0.3), color='black', alpha=0.8, size=0.5) + 
      geom_segment(data=d6[d6$strand=='+', ], aes(x=start/1e6, y=-1.3*grp, xend=end/1e6, yend=-1.3*grp), 
                  arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.5) +
      geom_segment(data=d6[d6$strand=='-', ], aes(x=end/1e6, y=-1.3*grp, xend=start/1e6, yend=-1.3*grp), 
                  arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.5) +
      geom_text_repel(data=d6, aes(x=(start+end)/2e6, y=-1.3*grp, label=gene_name), size=1.8, 
           box.padding = unit(0.15, "lines"), point.padding = unit(0.15, "lines")) + 
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        plot.title = element_text(size = 8),
                        axis.text = element_text(size=8,colour = "black"),
                        axis.title = element_text(size=8),
                        axis.title.x = element_text(vjust=0.2),
                        axis.title.y = element_text(vjust=0.2),
                        axis.ticks.length = unit(0.25,"cm"))
#dev.off()
# png(file=paste0('summary/plots/val_bid',d1$bid[1],'.png'),height=3,width=4,unit='in', res=300)
y1max <- ceiling(max(-log10(d3$LRT)))+1
y2max <- ceiling(max(c(d3$AFconc,d3$AFcase)))
p2 <- ggplot(data=d3, aes(x=hcut, y=-log10(LRT), color='LRT test')) +
     geom_line(size=0.6) + geom_point(size=0.3) +
     scale_y_continuous(name='-Log10(P)', limits = c(0, y1max+1), sec.axis=sec_axis(~.*(y2max/y1max),name='Frq per 10,000'))+
     scale_x_continuous(name='HDS error c', limits = c(0, 0.5), labels=seq(0,0.5,0.1)) +
     ggtitle('Validation in the DRIVE Imputed Data')+
     geom_line(data=d3, aes(x=hcut, y=AFcase*y1max/y2max, color='BCa'), size=0.6) +
     geom_line(data=d3, aes(x=hcut, y=AFconc*y1max/y2max, color='Control'), size=0.6) +
     scale_colour_manual(name='', breaks=c('LRT test','BCa','Control'), 
           values=c('LRT test'='black', 'BCa'="blue",'Control'="orange"))+
     geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black",alpha=0.5,size=0.4)+
     theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.title = element_text(size = 8),
                        axis.line = element_line(colour = "black"),
                        axis.text = element_text(size=8,colour = "black"),
                        axis.title = element_text(size=8),
                        axis.ticks.length = unit(0.25,"cm"),
                        legend.position=c(0.50,1), legend.direction='horizontal',
                        legend.key=element_blank(), legend.spacing.x = unit(0.0, 'cm'),
                        legend.text=element_text(size=5)) 
#dev.off()
p3 <- ggdraw()+draw_image(paste0('summary/plots/r2_ukb_bid',d1$bid,'.png'),scale=1)
p4 <- ggdraw()+draw_image(paste0('summary/plots/dpri_ukb_bid',d1$bid,'.png'),scale=1)
bottom_row <- plot_grid(p3, p4, p2, labels=c('B','C','D'), nrow=1, rel_widths=c(1,1,2), label_size=10)
pp[[idx]] <- plot_grid(p1, bottom_row, labels =c('A', ''), ncol=1, rel_heights=c(1.2,1),label_size=10)
}


pdf(file='summary/plots/assoc_plots.pdf', height=4,width=5)
pp[[1]]
pp[[2]]
pp[[3]]
pp[[4]]
pp[[5]]
pp[[6]]
pp[[7]]
pp[[8]]
pp[[9]]
pp[[10]]
pp[[11]]
dev.off()


