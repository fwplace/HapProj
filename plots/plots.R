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
