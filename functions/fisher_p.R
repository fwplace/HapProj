fisher_p <- function(z, y, mod){
  y <- as.numeric(y); z <- as.numeric(z)
  idx <- which(!is.na(y)&!is.na(z))
  yy <- y[idx]; zz <- z[idx]

  if(mod=='add'){
         c1 <- c(sum(yy==1&zz==0), sum(yy==1&zz==1), sum(yy==1&zz==2))
         c0 <- c(sum(yy==0&zz==0), sum(yy==0&zz==1), sum(yy==0&zz==2))
  }
  if(mod=='dom'){
         c1 <- c(sum(yy==1&zz==0), sum(yy==1&zz>=1))
         c0 <- c(sum(yy==0&zz==0), sum(yy==0&zz>=1))
  }
  if(mod=='rec'){
         c1 <- c(sum(yy==1&zz<=1), sum(yy==1&zz==2))
         c0 <- c(sum(yy==0&zz<=1), sum(yy==0&zz==2))
  }

 res <- data.frame(c0=paste(c0, collapse='/'), c1=paste(c1, collapse='/'),
                  p=tryCatch(fisher.test(rbind(c1, c0))$p.value, 
                  error=function(e){return(chisq.test(rbind(c1, c0))$p.value)}))

return(res)
}

