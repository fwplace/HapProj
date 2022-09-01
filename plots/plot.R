### Local association plots
library(ggplot2)
library(cowplot)
library(ggrepel)

load('../analysis/DRIVE4/rare13/tracks/snp.gwas.rdata')
load('../analysis/DRIVE4/rare13/tracks/recomb.female.rdata')
load('../analysis/DRIVE4/rare13/tracks/tad.HMEC.rdata')
load('../analysis/DRIVE4/rare13/tracks/loops.breast.rdata')
rare13.genes <- read.table(file='../analysis/DRIVE4/rare13/tracks/rare13.genes.v2.txt',header=T)
load('../analysis/DRIVE4/rare13/rare13.obj.rdata')
load('../analysis/DRIVE4/rare13/anno/vinfo.rdata')
d1 <- do.call('rbind', lapply(rare13.obj, function(x){x$sum}))
d2 <- do.call('rbind', lapply(rare13.obj, function(x){x$ukb}))
d1 <- d1[match(vinfo$hid, d1$hid), ]
d2 <- d2[match(vinfo$hid, d2$hid), ]
colnames(vinfo) <- gsub('comb.', 'ukb.', colnames(vinfo))
vinfo <- data.frame(d1[, c(1,2,4,11,7:9,13:18)], d2[, c(3:6,9:12)], vinfo[, -c(1:2)])
hdat <- vinfo[!duplicated(vinfo$bid), ]
anno <- read.csv(file='../analysis/DRIVE4/rare13/rare13.sum.csv')[, -1]

## bid
cex0=7
assoc_plots <- list()
dist0 <- 4e5
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
  labs <- paste0(d0$chr[1],':',min(d0$pos),'-',max(d0$pos))
  
  pp <- ggplot(data=d1, aes(x=pos/1e6, y=-log10(p))) + 
    geom_point(size=0.5, colour='black') +
    scale_x_continuous(name=paste0('Chr',hdat$chr[i], ', Mb'), limits = c(st/1e6, ed/1e6))+
    scale_y_continuous(name='-Log10(P)', limits = c(-16, ymax), breaks=seq(0,ymax,4), 
                       sec.axis=sec_axis(~.*(100/ymax), breaks=seq(0,100,20), name='cM/Mb'))+
    geom_line(data=d2, aes(x=pos/1e6, y=rate*ymax/100), color="brown", size=0.5, alpha=0.6)+
    ggtitle(paste0('Locus ',hdat$bid[i],': ', labs)) +
    geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "black", alpha=0.5,size=0.5)+
    geom_segment(data=data.frame(x1=min(d0$pos), x2=max(d0$pos), y=-log10(hdat$comb.cox[i])), 
                 aes(x=x1/1e6 , y = y, xend = x2/1e6, yend = y), color='orange', size=2) +
    geom_segment(data=d0, aes(x=pos/1e6, y = -log10(comb.cox)-0.3, xend = pos/1e6, 
                              yend = -log10(comb.cox)+0.3), color='black', size=0.5) +
    geom_segment(data=data.frame(x1=min(d0$pos[d0$validation==1]), 
                                 x2=max(d0$pos[d0$validation==1]), 
                                 y=-log10(hdat$p[i])), 
                 aes(x=x1/1e6 , y = y, xend = x2/1e6, yend = y), color='purple', size=2)+
    geom_segment(data=d0[d0$validation==1, ], 
                 aes(x=pos/1e6, y = -log10(p)-0.3, xend = pos/1e6, yend = -log10(p)+0.3), 
                 color='black', size=0.5)
  
  if(nrow(d3)>0) pp <- pp + geom_segment(data=d3, 
                                         aes(x=start/1e6 , y = -2, xend = end/1e6, yend = -2), 
                                         color='darkgreen', alpha=0.8, size=2)
  if(nrow(d4)>0) pp <- pp + geom_segment(data=d4, 
                                         aes(x=st1/1e6 , y = -4, xend = ed1/1e6, yend = -4), 
                                         color='magenta', size=2)+
    geom_segment(data=d4, 
                 aes(x=st2/1e6 , y = -4, xend = ed2/1e6, yend = -4), color='magenta', size=2)+
    geom_curve(data=d4, 
               aes(x=(st1+ed1)/2e6 , y = -4.1, xend = (st2+ed2)/2e6, yend = -4.1), 
               color='magenta', curvature = 0.4, size=0.3)
  
  if(sum(d5$strand=='+')>0) pp <- pp + 
    geom_segment(data=d5[d5$strand=='+', ],
                 aes(x=start/1e6, y=-3*group-7, xend=end/1e6, yend=-3*group-7), 
                 arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.5)
  
  if(sum(d5$strand=='-')>0) pp <- pp + 
    geom_segment(data=d5[d5$strand=='-', ], 
                 aes(x=end/1e6, y=-3*group-7, xend=start/1e6, yend=-3*group-7), 
                 arrow = arrow(length = unit(0.1, "cm")), color='blue', size=0.5)
  
  if(nrow(d5)>0) pp <- pp + 
    geom_text_repel(data=d5, aes(x=(start+end)/2e6, y=-3*group-7, label=gene_name), 
                    box.padding = 0.1, point.padding = 0.1, size=1.8)
  
  assoc_plots[[i]] <- pp + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black", size=0.5), 
          plot.title = element_text(face="bold", size = cex0+1),
          axis.text = element_text(size=cex0, colour = "black"), 
          axis.title = element_text(size=cex0),
          axis.title.x = element_text(vjust=-0.3),
          axis.title.y = element_text(vjust=-0.3),
          axis.ticks.length = unit(0.1,"cm"), 
          plot.margin = unit(c(0.1,0,0.1,0), "cm"))
}

df <- data.frame(x=rep(1:13,5), name=rep(anno$bid, 5),
                 type=rep(c('eQTL','Enhancer','TADBoundary','Chromatin Loop','Motif Change'), each=13), 
                 frq=c(anno$gtex8.bm/anno$w, anno$enhancer.breast/anno$w, anno$tad.boundary.20k/anno$w, 
                       anno$loop/anno$w, anno$motif.diff3/anno$w))
df$type <- factor(df$type, levels=c('eQTL','Enhancer','TADBoundary','Chromatin Loop','Motif Change'))
pp <- ggplot(data=df, aes(x=x, y=frq, col=type)) + 
  geom_line(size=0.6) +
  scale_x_continuous(name='Haplotype Locus', breaks=seq(1,13,1), labels = anno$bid)+
  scale_y_continuous(name='% Annotated Variants', limits=c(0,1), breaks=seq(0,1,0.2))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), panel.background = element_blank(),
                     axis.line = element_line(colour = "black"), 
                     plot.title = element_text(face="bold", size = cex0),
                     axis.text = element_text(size=cex0, colour = "black"), 
                     axis.title = element_text(size=cex0),
                     axis.title.x = element_text(vjust=0.6),axis.title.y = element_text(vjust=0.6),
                     axis.ticks.length = unit(0.25,"cm"), legend.title = element_blank()) 

pdf(file='Figure 2.pdf', height=9.5, width=7)
top <- plot_grid(assoc_plots[[1]], assoc_plots[[2]], assoc_plots[[3]],
                 assoc_plots[[4]],assoc_plots[[5]],assoc_plots[[6]],
                 assoc_plots[[7]],assoc_plots[[8]],assoc_plots[[9]], 
                 assoc_plots[[10]],assoc_plots[[11]],assoc_plots[[12]],
                 ncol=3, labels=LETTERS[1:12],label_size = 10, align = 'hv')
bottom <- plot_grid(assoc_plots[[13]], pp, rel_widths = c(1,2), 
                    nrow=1, labels=LETTERS[13:14], label_size = 10, align = 'hv')
plot_grid(top, bottom,rel_heights = c(4,1), ncol=1)
dev.off()




load('../analysis/DRIVE4/rare13/rare13.obj.rdata')
load('../analysis/DRIVE4/rare13/partition/rare13_ukbb_bigLD.rdata')
library(reshape)

ld_plots <- list()
for(idx in 1:13){
  dat <- rare13.obj[[idx]]$snp
  m1 <- as.matrix(rare13.obj[[idx]]$ld.ukbb)
  colnames(m1) <- rownames(m1) <- 1:nrow(dat)
  m2 <- melt(as.matrix(m1))
  res <- rare13.obj[[idx]]$sum
  labs <- paste0(dat$chr[1],':',min(dat$pos),'-',max(dat$pos))
  ld_plots[[idx]] <- ggplot(m2, aes(X1, X2)) + geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = "white", high = "red")+
    ggtitle(paste0('Locus ',res$bid[1],': ', labs)) +
    guides(fill=guide_legend(title=expression(D^"'"/~italic(r)^2)))+
    geom_segment(x=1, y = 1, xend = nrow(dat), yend = nrow(dat), color='black', alpha=0.5, size=0.6)+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.border = element_blank(), panel.background = element_blank(),
                       axis.line = element_blank(),
                       plot.title = element_text(size = 8, face='bold',margin=margin(0,0,0,0)),
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
    scale_x_continuous(name='Location of Variant', limits = c(0, nrow(dat)))+
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

pdf(file='Figure 3.pdf', height=8, width=7)
legend_pp <- get_legend(ld_plots[[1]])
plot_grid(ld_bigLD[[1]], ld_bigLD[[2]],ld_bigLD[[3]], ld_bigLD[[4]], ld_bigLD[[5]],ld_bigLD[[6]],
          ld_bigLD[[7]], ld_bigLD[[8]],ld_bigLD[[9]],ld_bigLD[[10]], ld_bigLD[[11]],ld_bigLD[[12]],
          ld_bigLD[[13]], legend_pp, ncol=3, labels=c(LETTERS[1:13],''), label_size = 10)
dev.off()



