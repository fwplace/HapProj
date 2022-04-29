logistic_p <- function(z, nz, y, x=NULL, get.lrt=FALSE, mod0.deviance=NA){
    if(nz==1) z <- as.numeric(z)
     y <- as.numeric(y)
    if(get.lrt){
        res <- as.data.frame(t(rep(NA,7)))
        colnames(res) <- c('beta','or', 'se', 'z', 'p','chi','p.lrt')
    }else{
        res <- as.data.frame(t(rep(NA,5)))
        colnames(res) <- c('beta','or', 'se', 'z', 'p')
    }

    dat <- cbind(z, y, x)
    is.error=FALSE
    fit <- tryCatch(glm(as.factor(y)~., data=as.data.frame(dat), family=binomial(link='logit')), error=function(e){is.error=TRUE})

    if(!is.error){
        if(get.lrt){
          if(is.na(mod0.deviance)){
               dat2 <- cbind(y,x)
               mod0.deviance <- glm(as.factor(y)~., data=as.data.frame(dat2), family=binomial(link='logit'))$deviance
          }
          chi <- mod0.deviance-fit$deviance
          p.lrt <- 1-pchisq(chi, df=nz)
        }
        if(nz==1){
               res <- as.data.frame(summary(fit)$coef)['z', ]
        }else{
               res <- as.data.frame(summary(fit)$coef)[colnames(z), ]
        }
        res <- cbind(res, exp(res[,1]))
        colnames(res) <- c('beta','se', 'z', 'p', 'or')
        res <- res[, c('beta','or', 'se', 'z', 'p')]
        if(get.lrt){
            res <- cbind(res, chi, p.lrt)
            colnames(res) <- c('beta','or', 'se', 'z', 'p','chi','p.lrt')
        }
    }
return(res)
}

