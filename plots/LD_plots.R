## LD plots
## intsall.packages('gaston')
setwd("~/fwork/ProjHap/breast/analysis/rare13")
library(gaston)
library(reshape)
library(ggplot2)
load('rare13.obj.rdata')
ld_plots <- list()
for(idx in 1:13){
  dat <- rare13.obj[[idx]]$snpset
  tmp <- which(dat$validation==1)
  m1 <- as.matrix(rare13.obj[[idx]]$ukb.ld[tmp, tmp])
  colnames(m1) <- rownames(m1) <- dat$snp[tmp]
  m2 <- melt(as.matrix(m1))
  res <- rare13.obj[[idx]]$sum
  
  if(idx==13){
     ld_plots[[idx]] <- ggplot(m2, aes(X1, X2)) + geom_tile(aes(fill = value)) +
      scale_fill_gradient(low = "white", high = "red")+
      ggtitle(paste0('Bid=',res$bid[1],': ', res$block[1])) +
      guides(fill=guide_legend(title=expression(D^"'"/~italic(r)^2)))+
       geom_abline(intercept = 0, slope = 1)+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.border = element_blank(), panel.background = element_blank(),
                         axis.line = element_blank(),plot.title = element_text(size = 10),
                         axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                         axis.text.y = element_text(size=10,colour = "black"), 
                         axis.title = element_blank())
  }else{
    ld_plots[[idx]] <- ggplot(m2, aes(X1, X2)) + geom_tile(aes(fill = value)) +
      scale_fill_gradient(low = "white", high = "red")+
      ggtitle(paste0('Bid=',res$bid[1],': ', res$block[1])) +
      guides(fill=guide_legend(title=expression(D^"'"/~italic(r)^2)))+
      geom_abline(intercept = 0, slope = 1)+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.border = element_blank(), panel.background = element_blank(),
                         axis.line = element_blank(),plot.title = element_text(size = 10),
                         axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                         axis.text.y = element_text(size=10,colour = "black"), 
                         axis.title = element_blank(),legend.position = 'none')
  }
  
}


pdf(file='Rare13_reduction_LD_plots.pdf',height=10,width=9)
legend_pp <- get_legend(ld_plots[[13]])
ld_plots[[13]] <- ld_plots[[13]] + theme(legend.position ='none')
plot_grid(ld_plots[[1]], ld_plots[[2]],ld_plots[[3]],
          ld_plots[[4]], ld_plots[[5]],ld_plots[[6]],
          ld_plots[[7]], ld_plots[[8]],ld_plots[[9]],
          ld_plots[[10]],ld_plots[[11]],ld_plots[[12]],
          ld_plots[[13]], legend_pp,ncol=3, labels=LETTERS[1:13])
dev.off()
