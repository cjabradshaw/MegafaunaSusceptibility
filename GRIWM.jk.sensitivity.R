#######################################################################################
## Gaussian-Resampled Inverse-Weighted McInerny (GRIWM) extinction estimator
## The model takes a .csv file with 2 columns as input headed [,1]=age and [,2]=sd
## GRIWM returns a variable named 'out' as an output such as:
## out[1] and out[3] = CI as quartiles at 97.5% and 2.5% of 10,000 iterations each iteration,
## being a resampling of the ages into their standard deviation
## out[2] = the median value = your final estimates
##
## Frédérik Saltré & Corey J. A. Bradshaw
## Flinders University
## March 2021
######################################################################################

rm(list=ls(all=TRUE))

GRIWM <- function(dat = dat, alpha = alpha, iter = iter) 
{
  itdiv <- iter/(iter/10);out<-matrix(0,1,3)
  dat <- dat[order(dat[,1],decreasing=F),1:2] # order data

  date4 <- dat[,1]; sd.vec <- dat[,2]
  k <- length(date4)
  T.up.vec <- T.mci.vec <- w.T.mci.vec <- rep(0,iter)
  T.up.vec <- T.mci.vec <- w.T.mci.vec <- rep(0,iter)
    
  for (c in 1:iter) {
    date.samp <- rep(0,k)
      
    for (b in 1:k) {
      date.samp[b] <- round(rnorm(1,date4[b],sd.vec[b]))}
      
    date.samp <- (sort(date.samp))

    ## calculate weighted McInerny date & confidence interval
    last.diff <- 1/(date.samp-date.samp[1])[-1]
    weight <- last.diff/last.diff[1]
      
    if (last.diff[1] == Inf) {
      weight <- last.diff/last.diff[2]
      weight <- weight[-1]}
      
    ldate <- length(date.samp)
    T.mci.lst.vec <- rep(0,ldate-1)
      
    for (m in 1:(ldate-1)) {
      date.it <- date.samp[1:(1+m)]
      date.age.it <- date.samp[1:(1+m)]
        
      date.mci.it <- rev(max(date.it) + 1 - date.it)
        
      k <- length(date.it)
      t.n <- date.mci.it[k]
      n <- k
      T.rng <- t.n - date.mci.it[1]
        
      i =t.n - t.n*log(alpha)/n # original version of GRIWM estimates
 
      T.mci.lst.vec[m] <- max(date.it) + 1 - i
    } 
            
    if (last.diff[1] == Inf) {
      w.T.mci.vec[c] <- round((sum(weight*T.mci.lst.vec[-1]))/sum(weight),0)}
      
    if (last.diff[1] != Inf) {
      w.T.mci.vec[c] <- round((sum(weight*T.mci.lst.vec))/sum(weight),0)}
      
    }
  T.wmci.vec.CI <- quantile(na.omit(w.T.mci.vec),probs = c(0.025, 0.5, 0.975))
  out[1]<-round(T.wmci.vec.CI[1], 0);out[2]<-round(T.wmci.vec.CI[2], 0);out[3]<-round(T.wmci.vec.CI[3], 0)
  return(out)
  rm(list=ls(all=TRUE))
} 

#########################################################################################################################################
# jack-knife sensitivity analysis
#########################################################################################################################################
## each 'XXXX.csv' is the two-column file of raw date estimates of the species' (XXXX) chronology

library(MASS) ## call libraries
dat <- read.table("XXXX.csv", header = T, sep = ",")
repet = 1000
alpha <- 0.05
mat <- matrix(0,repet,3);nbdat<-length(dat[,1])-5

for (i in 1:repet) {
  print(i)
  pkdel <- round(runif(1,0,nbdat))
  pkkeep <- length(dat[,1]) - pkdel
  sdatindex <- sort(sample(1:length(dat[,1]), pkkeep,replace = FALSE), decreasing = FALSE)
  dat2 <- dat[sdatindex,1:2]
  mat[i,] <- GRIWM(dat2, alpha, 10000)
  rm(dat2,pkdel,pkkeep,sdatindex)
}

head(mat)

