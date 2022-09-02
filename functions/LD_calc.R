LD_calc <- function(mat){ ## 2N x 2
   p1 <- mean(mat[,1])
   q1 <- mean(mat[,2])
   d <- mean(mat[,1]*mat[,2]) - p1*q1
   if(d<0)  dmax <- max(-1*p1*q1, -1*(1-p1)*(1-q1))
   if(d>=0) dmax <- min(p1*(1-q1), (1-p1)*q1)
   if(dmax==0){
         dpri <- r2 <- 0
   }else{
         dpri <- d/dmax
         r2  <- d^2/(p1*(1-p1)*q1*(1-q1))
   }
return(c(dpri,r2))
}