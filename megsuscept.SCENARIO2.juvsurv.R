##########################################################################################################################################
## megafauna demographic susceptibility model
## AIM: construct plausible stochastic demographic models for main Sahul megafauna to determine relative demographic
##      susceptibility to environmental change & novel predation (human) sources
##
##          VOMBATIFORM HERBIVORES: ✓Diprotodon (†), ✓Palorchestes (†), ✓Zygomaturus (†), ✓Phascolonus (†), ✓Vombatus ursinus
##          MACROPODIFORM HERBIVORES: ✓Protemnodon (†), ✓Osphranter rufus, ✓Sthenurus (†), ✓Simosthenurus (†), ✓Procoptodon (†), ✓Metasthenurus (†), ✓Notamacropus
##          LARGE BIRDS: ✓Genyornis (†), ✓Dromaius novaehollandiae, ✓Alectura lathami
##          CARNIVORES: ✓Sarcophilus, ✓Thylacinus (†), ✓Thylacoleo (†), ✓Dasyurus
##          MONOTREMES: ✓Megalibgwilia (†), ✓Tachyglossus
##
## Corey Bradshaw
## corey.bradshaw@flinders.edu.au
## Flinders University, August 2020
##########################################################################################################################################

###################################################################################
### *** run 'Sahul megafauna demographic susceptibility-base models.R' first *** ##
###################################################################################

#####################################################################
## SCENARIO 2: PROGRESSIVELY INCREASE JUVENILE (0-alpha) MORTALITY ##
#####################################################################

iter <- 10000
itdiv <- iter/10
Q.ext.thresh <- 100/2 # quasi-extinction threshold

juv.mort.inc.vec <- seq(0, 0.975, 0.025)

DP.Q.ext.pr <- PA.Q.ext.pr <- ZT.Q.ext.pr <-  PH.Q.ext.pr <- VU.Q.ext.pr <- PG.Q.ext.pr <- SS.Q.ext.pr <- PT.Q.ext.pr <- SO.Q.ext.pr <- MN.Q.ext.pr <- OR.Q.ext.pr <- NR.Q.ext.pr <- GN.Q.ext.pr <- DN.Q.ext.pr <- AL.Q.ext.pr <- TC.Q.ext.pr <- TH.Q.ext.pr <- SH.Q.ext.pr <- DM.Q.ext.pr <- TA.Q.ext.pr <- MR.Q.ext.pr <- rep(NA,length(juv.mort.inc.vec))

for (k in 1:length(juv.mort.inc.vec)) {

  ## DIPROTODON (DP)
  ## set storage matrices & vectors
  DP.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    DP.popmat <- DP.popmat.orig
    
    DP.n.mat <- matrix(0, nrow=DP.age.max+1,ncol=(t+1))
    DP.n.mat[,1] <- DP.init.vec

    for (i in 1:t) {
      # stochastic survival values
      DP.s.alpha <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$alpha
      DP.s.beta <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$beta
      DP.s.stoch <- rbeta(length(DP.s.alpha), DP.s.alpha, DP.s.beta)
      DP.s.stoch[1:DP.alpha] <- DP.s.stoch[1:DP.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/DP.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        DP.s.stoch <- DP.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
  
      # stochastic fertilty sampler (gaussian)
      DP.fert.stch <- rnorm(length(DP.popmat[,1]), DP.pred.p.mm, DP.m.sd.vec)

      DP.totN.i <- sum(DP.n.mat[,i], na.rm=T)
      DP.pred.red <- DP.a.lp/(1+(DP.totN.i/DP.b.lp)^DP.c.lp)
      
      diag(DP.popmat[2:(DP.age.max+1),]) <- (DP.s.stoch[-(DP.age.max+1)])*DP.pred.red
      DP.popmat[DP.age.max+1,DP.age.max+1] <- (DP.s.stoch[DP.age.max+1])*DP.pred.red
      DP.popmat[1,] <- ifelse(DP.fert.stch < 0, 0, DP.fert.stch)
      DP.n.mat[,i+1] <- DP.popmat %*% DP.n.mat[,i]

    } # end i loop
    
    DP.n.sums.mat[e,] <- (as.vector(colSums(DP.n.mat)))

  } # end e loop
  
  # total N
  DP.n.md <- apply(DP.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  DP.n.up <- apply(DP.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  DP.n.lo <- apply(DP.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  DP.Q.ext.mat <- ifelse(DP.n.sums.mat < Q.ext.thresh, 1, 0)
  DP.Q.ext.sum <- apply(DP.Q.ext.mat[,ceiling(DP.gen.l):dim(DP.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  DP.Q.ext.pr[k] <- length(which(DP.Q.ext.sum > 0)) / iter
  print("DIPROTODON")
  

    
  ## PALORCHESTES (PA)
  ## set storage matrices & vectors
  PA.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    PA.popmat <- PA.popmat.orig
    
    PA.n.mat <- matrix(0, nrow=PA.age.max+1,ncol=(t+1))
    PA.n.mat[,1] <- PA.init.vec

    for (i in 1:t) {
      # stochastic survival values
      PA.s.alpha <- estBetaParams(PA.Sx, PA.s.sd.vec^2)$alpha
      PA.s.beta <- estBetaParams(PA.Sx, PA.s.sd.vec^2)$beta
      PA.s.stoch <- rbeta(length(PA.s.alpha), PA.s.alpha, PA.s.beta)
      PA.s.stoch[1:PA.alpha] <- PA.s.stoch[1:PA.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/PA.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PA.s.stoch <- PA.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      PA.fert.stch <- rnorm(length(PA.popmat[,1]), PA.pred.p.mm, PA.m.sd.vec)
      
      PA.totN.i <- sum(PA.n.mat[,i], na.rm=T)
      PA.pred.red <- PA.a.lp/(1+(PA.totN.i/PA.b.lp)^PA.c.lp)
      
      diag(PA.popmat[2:(PA.age.max+1),]) <- (PA.s.stoch[-(PA.age.max+1)])*PA.pred.red
      PA.popmat[PA.age.max+1,PA.age.max+1] <- (PA.s.stoch[PA.age.max+1])*PA.pred.red
      PA.popmat[1,] <- ifelse(PA.fert.stch < 0, 0, PA.fert.stch)
      PA.n.mat[,i+1] <- PA.popmat %*% PA.n.mat[,i]
      
    } # end i loop
    
    PA.n.sums.mat[e,] <- (as.vector(colSums(PA.n.mat)))

  } # end e loop
  
  # total N
  PA.n.md <- apply(PA.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PA.n.up <- apply(PA.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PA.n.lo <- apply(PA.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  PA.Q.ext.mat <- ifelse(PA.n.sums.mat < Q.ext.thresh, 1, 0)
  PA.Q.ext.sum <- apply(PA.Q.ext.mat[,ceiling(PA.gen.l):dim(PA.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PA.Q.ext.pr[k] <- length(which(PA.Q.ext.sum > 0)) / iter
  print("PALORCHESTES")
  

    
  ## ZYGOMATURUS (ZT)
  ## set storage matrices & vectors
  ZT.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    ZT.popmat <- ZT.popmat.orig
    
    ZT.n.mat <- matrix(0, nrow=ZT.age.max+1,ncol=(t+1))
    ZT.n.mat[,1] <- ZT.init.vec

    for (i in 1:t) {
      # stochastic survival values
      ZT.s.alpha <- estBetaParams(ZT.Sx, ZT.s.sd.vec^2)$alpha
      ZT.s.beta <- estBetaParams(ZT.Sx, ZT.s.sd.vec^2)$beta
      ZT.s.stoch <- rbeta(length(ZT.s.alpha), ZT.s.alpha, ZT.s.beta)
      ZT.s.stoch[1:ZT.alpha] <- ZT.s.stoch[1:ZT.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/ZT.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        ZT.s.stoch <- ZT.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      ZT.fert.stch <- rnorm(length(ZT.popmat[,1]), ZT.pred.p.mm, ZT.m.sd.vec)
      
      ZT.totN.i <- sum(ZT.n.mat[,i], na.rm=T)
      ZT.pred.red <- ZT.a.lp/(1+(ZT.totN.i/ZT.b.lp)^ZT.c.lp)
      
      diag(ZT.popmat[2:(ZT.age.max+1),]) <- (ZT.s.stoch[-(ZT.age.max+1)])*ZT.pred.red
      ZT.popmat[ZT.age.max+1,ZT.age.max+1] <- (ZT.s.stoch[ZT.age.max+1])*ZT.pred.red
      ZT.popmat[1,] <- ifelse(ZT.fert.stch < 0, 0, ZT.fert.stch)
      ZT.n.mat[,i+1] <- ZT.popmat %*% ZT.n.mat[,i]

    } # end i loop
    
    ZT.n.sums.mat[e,] <- (as.vector(colSums(ZT.n.mat)))

  } # end e loop
  
  # total N
  ZT.n.md <- apply(ZT.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  ZT.n.up <- apply(ZT.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  ZT.n.lo <- apply(ZT.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  ZT.Q.ext.mat <- ifelse(ZT.n.sums.mat < Q.ext.thresh, 1, 0)
  ZT.Q.ext.sum <- apply(ZT.Q.ext.mat[,ceiling(ZT.gen.l):dim(ZT.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  ZT.Q.ext.pr[k] <- length(which(ZT.Q.ext.sum > 0)) / iter
  print("ZYGOMATURUS")
  

  
  ## PHASCOLONUS (PH)
  ## set storage matrices & vectors
  PH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    PH.popmat <- PH.popmat.orig
    
    PH.n.mat <- matrix(0, nrow=PH.age.max+1,ncol=(t+1))
    PH.n.mat[,1] <- PH.init.vec

    for (i in 1:t) {
      # stochastic survival values
      PH.s.alpha <- estBetaParams(PH.Sx, PH.s.sd.vec^2)$alpha
      PH.s.beta <- estBetaParams(PH.Sx, PH.s.sd.vec^2)$beta
      PH.s.stoch <- rbeta(length(PH.s.alpha), PH.s.alpha, PH.s.beta)
      PH.s.stoch[1:PH.alpha] <- PH.s.stoch[1:PH.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/PH.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PH.s.stoch <- PH.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      PH.fert.stch <- rnorm(length(PH.popmat[,1]), PH.pred.p.mm, PH.m.sd.vec)
      
      PH.totN.i <- sum(PH.n.mat[,i], na.rm=T)
      PH.pred.red <- PH.a.lp/(1+(PH.totN.i/PH.b.lp)^PH.c.lp)
      
      diag(PH.popmat[2:(PH.age.max+1),]) <- (PH.s.stoch[-(PH.age.max+1)])*PH.pred.red
      PH.popmat[PH.age.max+1,PH.age.max+1] <- (PH.s.stoch[PH.age.max+1])*PH.pred.red
      PH.popmat[1,] <- ifelse(PH.fert.stch < 0, 0, PH.fert.stch)
      PH.n.mat[,i+1] <- PH.popmat %*% PH.n.mat[,i]

    } # end i loop
    
    PH.n.sums.mat[e,] <- (as.vector(colSums(PH.n.mat)))

    } # end e loop
  
  # total N
  PH.n.md <- apply(PH.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PH.n.up <- apply(PH.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PH.n.lo <- apply(PH.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  PH.Q.ext.mat <- ifelse(PH.n.sums.mat < Q.ext.thresh, 1, 0)
  PH.Q.ext.sum <- apply(PH.Q.ext.mat[,ceiling(PH.gen.l):dim(PH.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PH.Q.ext.pr[k] <- length(which(PH.Q.ext.sum > 0)) / iter
  print("PHASCOLONUS")

  
  
  ## VOMBATUS (VU)
  ## set storage matrices & vectors
  VU.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    VU.popmat <- VU.popmat.orig
    
    VU.n.mat <- matrix(0, nrow=VU.age.max+1,ncol=(t+1))
    VU.n.mat[,1] <- VU.init.vec

    for (i in 1:t) {
      # stochastic survival values
      VU.s.alpha <- estBetaParams(VU.Sx, VU.s.sd.vec^2)$alpha
      VU.s.beta <- estBetaParams(VU.Sx, VU.s.sd.vec^2)$beta
      VU.s.stoch <- rbeta(length(VU.s.alpha), VU.s.alpha, VU.s.beta)
      VU.s.stoch[1:VU.alpha] <- VU.s.stoch[1:VU.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/VU.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        VU.s.stoch <- VU.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      VU.fert.stch <- rnorm(length(VU.popmat[,1]), VU.pred.p.mm, VU.m.sd.vec)
      
      VU.totN.i <- sum(VU.n.mat[,i], na.rm=T)
      VU.pred.red <- VU.a.lp/(1+(VU.totN.i/VU.b.lp)^VU.c.lp)
      
      diag(VU.popmat[2:(VU.age.max+1),]) <- (VU.s.stoch[-(VU.age.max+1)])*VU.pred.red
      VU.popmat[VU.age.max+1,VU.age.max+1] <- 0 # (VU.s.stoch[VU.age.max+1])*VU.pred.red
      VU.popmat[1,] <- ifelse(VU.fert.stch < 0, 0, VU.fert.stch)
      VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]

    } # end i loop
    
    VU.n.sums.mat[e,] <- (as.vector(colSums(VU.n.mat)))

  } # end e loop
  
  # total N
  VU.n.md <- apply(VU.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  VU.n.up <- apply(VU.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  VU.n.lo <- apply(VU.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  VU.Q.ext.mat <- ifelse(VU.n.sums.mat < Q.ext.thresh, 1, 0)
  VU.Q.ext.sum <- apply(VU.Q.ext.mat[,ceiling(VU.gen.l):dim(VU.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  VU.Q.ext.pr[k] <- length(which(VU.Q.ext.sum > 0)) / iter
  print("VOMBATUS")
  

  
  ## PROCOPTODON (PG)
  ## set storage matrices & vectors
  PG.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    PG.popmat <- PG.popmat.orig
    
    PG.n.mat <- matrix(0, nrow=PG.age.max+1,ncol=(t+1))
    PG.n.mat[,1] <- PG.init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      PG.s.alpha <- estBetaParams(PG.Sx, PG.s.sd.vec^2)$alpha
      PG.s.beta <- estBetaParams(PG.Sx, PG.s.sd.vec^2)$beta
      PG.s.stoch <- rbeta(length(PG.s.alpha), PG.s.alpha, PG.s.beta)
      PG.s.stoch[1:PG.alpha] <- PG.s.stoch[1:PG.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/PG.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PG.s.stoch <- PG.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      PG.fert.stch <- rnorm(length(PG.popmat[,1]), PG.pred.p.mm, PG.m.sd.vec)
      
      PG.totN.i <- sum(PG.n.mat[,i], na.rm=T)
      PG.pred.red <- PG.a.lp/(1+(PG.totN.i/PG.b.lp)^PG.c.lp)
      
      diag(PG.popmat[2:(PG.age.max+1),]) <- (PG.s.stoch[-(PG.age.max+1)])*PG.pred.red
      PG.popmat[PG.age.max+1,PG.age.max+1] <- (PG.s.stoch[PG.age.max+1])*PG.pred.red
      PG.popmat[1,] <- ifelse(PG.fert.stch < 0, 0, PG.fert.stch)
      PG.n.mat[,i+1] <- PG.popmat %*% PG.n.mat[,i]

    } # end i loop
    
    PG.n.sums.mat[e,] <- (as.vector(colSums(PG.n.mat)))

  } # end e loop
  
  # total N
  PG.n.md <- apply(PG.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PG.n.up <- apply(PG.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PG.n.lo <- apply(PG.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  PG.Q.ext.mat <- ifelse(PG.n.sums.mat < Q.ext.thresh, 1, 0)
  PG.Q.ext.sum <- apply(PG.Q.ext.mat[,ceiling(PG.gen.l):dim(PG.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PG.Q.ext.pr[k] <- length(which(PG.Q.ext.sum > 0)) / iter
  print("PROCOPTODON")
  
  

  ## STHENURUS (SS)
  ## set storage matrices & vectors
  SS.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    SS.popmat <- SS.popmat.orig
    
    SS.n.mat <- matrix(0, nrow=SS.age.max+1,ncol=(t+1))
    SS.n.mat[,1] <- SS.init.vec

    for (i in 1:t) {
      # stochastic survival values
      SS.s.alpha <- estBetaParams(SS.Sx, SS.s.sd.vec^2)$alpha
      SS.s.beta <- estBetaParams(SS.Sx, SS.s.sd.vec^2)$beta
      SS.s.stoch <- rbeta(length(SS.s.alpha), SS.s.alpha, SS.s.beta)
      SS.s.stoch[1:SS.alpha] <- SS.s.stoch[1:SS.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/SS.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        SS.s.stoch <- SS.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      SS.fert.stch <- rnorm(length(SS.popmat[,1]), SS.pred.p.mm, SS.m.sd.vec)
      
      SS.totN.i <- sum(SS.n.mat[,i], na.rm=T)
      SS.pred.red <- SS.a.lp/(1+(SS.totN.i/SS.b.lp)^SS.c.lp)
      
      diag(SS.popmat[2:(SS.age.max+1),]) <- (SS.s.stoch[-(SS.age.max+1)])*SS.pred.red
      SS.popmat[SS.age.max+1,SS.age.max+1] <- (SS.s.stoch[SS.age.max+1])*SS.pred.red
      SS.popmat[1,] <- ifelse(SS.fert.stch < 0, 0, SS.fert.stch)
      SS.n.mat[,i+1] <- SS.popmat %*% SS.n.mat[,i]

    } # end i loop
    
    SS.n.sums.mat[e,] <- (as.vector(colSums(SS.n.mat)))
 
  } # end e loop
  
  # total N
  SS.n.md <- apply(SS.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  SS.n.up <- apply(SS.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  SS.n.lo <- apply(SS.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  SS.Q.ext.mat <- ifelse(SS.n.sums.mat < Q.ext.thresh, 1, 0)
  SS.Q.ext.sum <- apply(SS.Q.ext.mat[,ceiling(SS.gen.l):dim(SS.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  SS.Q.ext.pr[k] <- length(which(SS.Q.ext.sum > 0)) / iter
  print("STHENURUS")
  
  
  
  ## PROTEMNODON (PT)
  ## set storage matrices & vectors
  PT.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    PT.popmat <- PT.popmat.orig
    
    PT.n.mat <- matrix(0, nrow=PT.age.max+1,ncol=(t+1))
    PT.n.mat[,1] <- PT.init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      PT.s.alpha <- estBetaParams(PT.Sx, PT.s.sd.vec^2)$alpha
      PT.s.beta <- estBetaParams(PT.Sx, PT.s.sd.vec^2)$beta
      PT.s.stoch <- rbeta(length(PT.s.alpha), PT.s.alpha, PT.s.beta)
      PT.s.stoch[1:PT.alpha] <- PT.s.stoch[1:PT.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/PT.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PT.s.stoch <- PT.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      PT.fert.stch <- rnorm(length(PT.popmat[,1]), PT.pred.p.mm, PT.m.sd.vec)
      
      PT.totN.i <- sum(PT.n.mat[,i], na.rm=T)
      PT.pred.red <- PT.a.lp/(1+(PT.totN.i/PT.b.lp)^PT.c.lp)
      
      diag(PT.popmat[2:(PT.age.max+1),]) <- (PT.s.stoch[-(PT.age.max+1)])*PT.pred.red
      PT.popmat[PT.age.max+1,PT.age.max+1] <- (PT.s.stoch[PT.age.max+1])*PT.pred.red
      PT.popmat[1,] <- ifelse(PT.fert.stch < 0, 0, PT.fert.stch)
      PT.n.mat[,i+1] <- PT.popmat %*% PT.n.mat[,i]
      
    } # end i loop
    
    PT.n.sums.mat[e,] <- (as.vector(colSums(PT.n.mat)))
    
  } # end e loop
  
  # total N
  PT.n.md <- apply(PT.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PT.n.up <- apply(PT.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PT.n.lo <- apply(PT.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  PT.Q.ext.mat <- ifelse(PT.n.sums.mat < Q.ext.thresh, 1, 0)
  PT.Q.ext.sum <- apply(PT.Q.ext.mat[,ceiling(PT.gen.l):dim(PT.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PT.Q.ext.pr[k] <- length(which(PT.Q.ext.sum > 0)) / iter
  print("PROTEMNODON")

  
  
  ## SIMOSTHENURUS (SO)
  ## set storage matrices & vectors
  SO.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    SO.popmat <- SO.popmat.orig
    
    SO.n.mat <- matrix(0, nrow=SO.age.max+1,ncol=(t+1))
    SO.n.mat[,1] <- SO.init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      SO.s.alpha <- estBetaParams(SO.Sx, SO.s.sd.vec^2)$alpha
      SO.s.beta <- estBetaParams(SO.Sx, SO.s.sd.vec^2)$beta
      SO.s.stoch <- rbeta(length(SO.s.alpha), SO.s.alpha, SO.s.beta)
      SO.s.stoch[1:SO.alpha] <- SO.s.stoch[1:SO.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/SO.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        SO.s.stoch <- SO.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      SO.fert.stch <- rnorm(length(SO.popmat[,1]), SO.pred.p.mm, SO.m.sd.vec)
      
      SO.totN.i <- sum(SO.n.mat[,i], na.rm=T)
      SO.pred.red <- SO.a.lp/(1+(SO.totN.i/SO.b.lp)^SO.c.lp)
      
      diag(SO.popmat[2:(SO.age.max+1),]) <- (SO.s.stoch[-(SO.age.max+1)])*SO.pred.red
      SO.popmat[SO.age.max+1,SO.age.max+1] <- (SO.s.stoch[SO.age.max+1])*SO.pred.red
      SO.popmat[1,] <- ifelse(SO.fert.stch < 0, 0, SO.fert.stch)
      SO.n.mat[,i+1] <- SO.popmat %*% SO.n.mat[,i]
      
    } # end i loop
    
    SO.n.sums.mat[e,] <- (as.vector(colSums(SO.n.mat)))
    
  } # end e loop
  
  # total N
  SO.n.md <- apply(SO.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  SO.n.up <- apply(SO.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  SO.n.lo <- apply(SO.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  SO.Q.ext.mat <- ifelse(SO.n.sums.mat < Q.ext.thresh, 1, 0)
  SO.Q.ext.sum <- apply(SO.Q.ext.mat[,ceiling(SO.gen.l):dim(SO.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  SO.Q.ext.pr[k] <- length(which(SO.Q.ext.sum > 0)) / iter
  print("SIMOSTHENURUS")

  
  
  ## METASTHENURUS (MN)
  ## set storage matrices & vectors
  MN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    MN.popmat <- MN.popmat.orig
    
    MN.n.mat <- matrix(0, nrow=MN.age.max+1,ncol=(t+1))
    MN.n.mat[,1] <- MN.init.vec

    for (i in 1:t) {
      # stochastic survival values
      MN.s.alpha <- estBetaParams(MN.Sx, MN.s.sd.vec^2)$alpha
      MN.s.beta <- estBetaParams(MN.Sx, MN.s.sd.vec^2)$beta
      MN.s.stoch <- rbeta(length(MN.s.alpha), MN.s.alpha, MN.s.beta)
      MN.s.stoch[1:MN.alpha] <- MN.s.stoch[1:MN.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/MN.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        MN.s.stoch <- MN.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      MN.fert.stch <- rnorm(length(MN.popmat[,1]), MN.pred.p.mm, MN.m.sd.vec)
      
      MN.totN.i <- sum(MN.n.mat[,i], na.rm=T)
      MN.pred.red <- MN.a.lp/(1+(MN.totN.i/MN.b.lp)^MN.c.lp)
      
      diag(MN.popmat[2:(MN.age.max+1),]) <- (MN.s.stoch[-(MN.age.max+1)])*MN.pred.red
      MN.popmat[MN.age.max+1,MN.age.max+1] <- (MN.s.stoch[MN.age.max+1])*MN.pred.red
      MN.popmat[1,] <- ifelse(MN.fert.stch < 0, 0, MN.fert.stch)
      MN.n.mat[,i+1] <- MN.popmat %*% MN.n.mat[,i]
      
    } # end i loop
    
    MN.n.sums.mat[e,] <- (as.vector(colSums(MN.n.mat)))

  } # end e loop
  
  # total N
  MN.n.md <- apply(MN.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  MN.n.up <- apply(MN.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  MN.n.lo <- apply(MN.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  MN.Q.ext.mat <- ifelse(MN.n.sums.mat < Q.ext.thresh, 1, 0)
  MN.Q.ext.sum <- apply(MN.Q.ext.mat[,ceiling(MN.gen.l):dim(MN.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  MN.Q.ext.pr[k] <- length(which(MN.Q.ext.sum > 0)) / iter
  print("METASTHENURUS")

  
  
  ## OSPHRANTER (OR)
  ## set storage matrices & vectors
  OR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    OR.popmat <- OR.popmat.orig
    
    OR.n.mat <- matrix(0, nrow=OR.age.max+1,ncol=(t+1))
    OR.n.mat[,1] <- OR.init.vec

    for (i in 1:t) {
      # stochastic survival values
      OR.s.alpha <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$alpha
      OR.s.beta <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$beta
      OR.s.stoch <- rbeta(length(OR.s.alpha), OR.s.alpha, OR.s.beta)
      OR.s.stoch[1:OR.alpha] <- OR.s.stoch[1:OR.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/OR.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        OR.s.stoch <- OR.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      OR.fert.stch <- rnorm(length(OR.popmat[,1]), OR.pred.p.mm, OR.m.sd.vec)
      
      OR.totN.i <- sum(OR.n.mat[,i], na.rm=T)
      OR.pred.red <- OR.a.lp/(1+(OR.totN.i/OR.b.lp)^OR.c.lp)
      
      diag(OR.popmat[2:(OR.age.max+1),]) <- (OR.s.stoch[-(OR.age.max+1)])*OR.pred.red
      OR.popmat[OR.age.max+1,OR.age.max+1] <- (OR.s.stoch[OR.age.max+1])*OR.pred.red
      OR.popmat[1,] <- ifelse(OR.fert.stch < 0, 0, OR.fert.stch)
      OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]

    } # end i loop
    
    OR.n.sums.mat[e,] <- (as.vector(colSums(OR.n.mat)))

  } # end e loop
  
  # total N
  OR.n.md <- apply(OR.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  OR.n.up <- apply(OR.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  OR.n.lo <- apply(OR.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  OR.Q.ext.mat <- ifelse(OR.n.sums.mat < Q.ext.thresh, 1, 0)
  OR.Q.ext.sum <- apply(OR.Q.ext.mat[,ceiling(OR.gen.l):dim(OR.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  OR.Q.ext.pr[k] <- length(which(OR.Q.ext.sum > 0)) / iter
  print("OSPHRANTER")
  
  
  
  ## NOTAMACROPUS rufogriseus (NR)
  ## set storage matrices & vectors
  NR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    NR.popmat <- NR.popmat.orig
    
    NR.n.mat <- matrix(0, nrow=NR.age.max+1,ncol=(t+1))
    NR.n.mat[,1] <- NR.init.vec

    for (i in 1:t) {
      # stochastic survival values
      NR.s.alpha <- estBetaParams(NR.Sx, NR.s.sd.vec^2)$alpha
      NR.s.beta <- estBetaParams(NR.Sx, NR.s.sd.vec^2)$beta
      NR.s.stoch <- rbeta(length(NR.s.alpha), NR.s.alpha, NR.s.beta)
      NR.s.stoch[1:NR.alpha] <- NR.s.stoch[1:NR.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/NR.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        NR.s.stoch <- NR.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      NR.fert.stch <- rnorm(length(NR.popmat[,1]), NR.pred.p.mm, NR.m.sd.vec)
      
      NR.totN.i <- sum(NR.n.mat[,i], na.rm=T)
      NR.pred.red <- NR.a.lp/(1+(NR.totN.i/NR.b.lp)^NR.c.lp)
      
      diag(NR.popmat[2:(NR.age.max+1),]) <- (NR.s.stoch[-(NR.age.max+1)])*NR.pred.red
      NR.popmat[NR.age.max+1,NR.age.max+1] <- (NR.s.stoch[NR.age.max+1])*NR.pred.red
      NR.popmat[1,] <- ifelse(NR.fert.stch < 0, 0, NR.fert.stch)
      NR.n.mat[,i+1] <- NR.popmat %*% NR.n.mat[,i]

    } # end i loop
    
    NR.n.sums.mat[e,] <- (as.vector(colSums(NR.n.mat)))

  } # end e loop
  
  # total N
  NR.n.md <- apply(NR.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  NR.n.up <- apply(NR.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  NR.n.lo <- apply(NR.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  NR.Q.ext.mat <- ifelse(NR.n.sums.mat < Q.ext.thresh, 1, 0)
  NR.Q.ext.sum <- apply(NR.Q.ext.mat[,ceiling(NR.gen.l):dim(NR.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  NR.Q.ext.pr[k] <- length(which(NR.Q.ext.sum > 0)) / iter
  print("NOTAMACROPUS")
  

    
  ## GENYORNIS (GN)
  ## set storage matrices & vectors
  GN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    GN.popmat <- GN.popmat.orig
    
    GN.n.mat <- matrix(0, nrow=GN.age.max+1,ncol=(t+1))
    GN.n.mat[,1] <- GN.init.vec

    for (i in 1:t) {
      # stochastic survival values
      GN.s.alpha <- estBetaParams(GN.Sx, GN.s.sd.vec^2)$alpha
      GN.s.beta <- estBetaParams(GN.Sx, GN.s.sd.vec^2)$beta
      GN.s.stoch <- rbeta(length(GN.s.alpha), GN.s.alpha, GN.s.beta)
      GN.s.stoch[1:GN.alpha] <- GN.s.stoch[1:GN.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/GN.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        GN.s.stoch <- GN.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      GN.fert.stch <- rnorm(length(GN.popmat[,1]), GN.pred.p.mm, GN.m.sd.vec)
      
      GN.totN.i <- sum(GN.n.mat[,i], na.rm=T)
      GN.pred.red <- GN.a.lp/(1+(GN.totN.i/GN.b.lp)^GN.c.lp)
      
      diag(GN.popmat[2:(GN.age.max+1),]) <- (GN.s.stoch[-(GN.age.max+1)])*GN.pred.red
      GN.popmat[GN.age.max+1,GN.age.max+1] <- (GN.s.stoch[GN.age.max+1])*GN.pred.red
      GN.popmat[1,] <- ifelse(GN.fert.stch < 0, 0, GN.fert.stch)
      GN.n.mat[,i+1] <- GN.popmat %*% GN.n.mat[,i]

    } # end i loop
    
    GN.n.sums.mat[e,] <- (as.vector(colSums(GN.n.mat)))

  } # end e loop
  
  # total N
  GN.n.md <- apply(GN.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  GN.n.up <- apply(GN.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  GN.n.lo <- apply(GN.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  GN.Q.ext.mat <- ifelse(GN.n.sums.mat < Q.ext.thresh, 1, 0)
  GN.Q.ext.sum <- apply(GN.Q.ext.mat[,ceiling(GN.gen.l):dim(GN.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  GN.Q.ext.pr[k] <- length(which(GN.Q.ext.sum > 0)) / iter
  print("GENYORNIS")
  
  
  
  ## DROMAIUS (GN)
  ## set storage matrices & vectors
  DN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    DN.popmat <- DN.popmat.orig
    
    DN.n.mat <- matrix(0, nrow=DN.age.max+1,ncol=(t+1))
    DN.n.mat[,1] <- DN.init.vec

    for (i in 1:t) {
      # stochastic survival values
      DN.s.alpha <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$alpha
      DN.s.beta <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$beta
      DN.s.stoch <- rbeta(length(DN.s.alpha), DN.s.alpha, DN.s.beta)
      DN.s.stoch[1:DN.alpha] <- DN.s.stoch[1:DN.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/DN.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        DN.s.stoch <- DN.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      DN.fert.stch <- rnorm(length(DN.popmat[,1]), DN.pred.p.mm, DN.m.sd.vec)
      
      DN.totN.i <- sum(DN.n.mat[,i], na.rm=T)
      DN.pred.red <- DN.a.lp/(1+(DN.totN.i/DN.b.lp)^DN.c.lp)
      
      diag(DN.popmat[2:(DN.age.max+1),]) <- (DN.s.stoch[-(DN.age.max+1)])*DN.pred.red
      DN.popmat[DN.age.max+1,DN.age.max+1] <- (DN.s.stoch[DN.age.max+1])*DN.pred.red
      DN.popmat[1,] <- ifelse(DN.fert.stch < 0, 0, DN.fert.stch)
      DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]

    } # end i loop
    
    DN.n.sums.mat[e,] <- (as.vector(colSums(DN.n.mat)))

  } # end e loop
  
  # total N
  DN.n.md <- apply(DN.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  DN.n.up <- apply(DN.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  DN.n.lo <- apply(DN.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  DN.Q.ext.mat <- ifelse(DN.n.sums.mat < Q.ext.thresh, 1, 0)
  DN.Q.ext.sum <- apply(DN.Q.ext.mat[,ceiling(DN.gen.l):dim(DN.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  DN.Q.ext.pr[k] <- length(which(DN.Q.ext.sum > 0)) / iter
  print("DROMAIUS")
  
  

  ## ALECTURA (AL)
  ## set storage matrices & vectors
  AL.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    AL.popmat <- AL.popmat.orig
    
    AL.n.mat <- matrix(0, nrow=AL.age.max+1,ncol=(t+1))
    AL.n.mat[,1] <- AL.init.vec

    for (i in 1:t) {
      # stochastic survival values
      AL.s.alpha <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$alpha
      AL.s.beta <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$beta
      AL.s.stoch <- rbeta(length(AL.s.alpha), AL.s.alpha, AL.s.beta)
      AL.s.stoch[1:AL.alpha] <- AL.s.stoch[1:AL.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/AL.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        AL.s.stoch <- AL.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      AL.fert.stch <- rnorm(length(AL.popmat[,1]), AL.pred.p.mm, AL.m.sd.vec)
      
      AL.totN.i <- sum(AL.n.mat[,i], na.rm=T)
      AL.pred.red <- AL.a.lp/(1+(AL.totN.i/AL.b.lp)^AL.c.lp)
      
      diag(AL.popmat[2:(AL.age.max+1),]) <- (AL.s.stoch[-(AL.age.max+1)])*AL.pred.red
      AL.popmat[AL.age.max+1,AL.age.max+1] <- (AL.s.stoch[AL.age.max+1])*AL.pred.red
      AL.popmat[1,] <- ifelse(AL.fert.stch < 0, 0, AL.fert.stch)
      AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]

    } # end i loop
    
    AL.n.sums.mat[e,] <- (as.vector(colSums(AL.n.mat)))

  } # end e loop
  
  # total N
  AL.n.md <- apply(AL.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  AL.n.up <- apply(AL.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  AL.n.lo <- apply(AL.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  AL.Q.ext.mat <- ifelse(AL.n.sums.mat < Q.ext.thresh, 1, 0)
  AL.Q.ext.sum <- apply(AL.Q.ext.mat[,ceiling(AL.gen.l):dim(AL.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  AL.Q.ext.pr[k] <- length(which(AL.Q.ext.sum > 0)) / iter
  print("ALECTURA")
  
  
  
  ## THYLACOLEO (TC)
  ## set storage matrices & vectors
  TC.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    TC.popmat <- TC.popmat.orig
    
    TC.n.mat <- matrix(0, nrow=TC.age.max+1,ncol=(t+1))
    TC.n.mat[,1] <- TC.init.vec

    for (i in 1:t) {
      # stochastic survival values
      TC.s.alpha <- estBetaParams(TC.Sx, TC.s.sd.vec^2)$alpha
      TC.s.beta <- estBetaParams(TC.Sx, TC.s.sd.vec^2)$beta
      TC.s.stoch <- rbeta(length(TC.s.alpha), TC.s.alpha, TC.s.beta)
      TC.s.stoch[1:TC.alpha] <- TC.s.stoch[1:TC.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/TC.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        TC.s.stoch <- TC.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      TC.fert.stch <- rnorm(length(TC.popmat[,1]), TC.pred.p.mm, TC.m.sd.vec)
      
      TC.totN.i <- sum(TC.n.mat[,i], na.rm=T)
      TC.pred.red <- TC.a.lp/(1+(TC.totN.i/TC.b.lp)^TC.c.lp)
      
      diag(TC.popmat[2:(TC.age.max+1),]) <- (TC.s.stoch[-(TC.age.max+1)])*TC.pred.red
      TC.popmat[TC.age.max+1,TC.age.max+1] <- (TC.s.stoch[TC.age.max+1])*TC.pred.red
      TC.popmat[1,] <- ifelse(TC.fert.stch < 0, 0, TC.fert.stch)
      TC.n.mat[,i+1] <- TC.popmat %*% TC.n.mat[,i]

    } # end i loop
    
    TC.n.sums.mat[e,] <- (as.vector(colSums(TC.n.mat)))

  } # end e loop
  
  # total N
  TC.n.md <- apply(TC.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  TC.n.up <- apply(TC.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  TC.n.lo <- apply(TC.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  TC.Q.ext.mat <- ifelse(TC.n.sums.mat < Q.ext.thresh, 1, 0)
  TC.Q.ext.sum <- apply(TC.Q.ext.mat[,ceiling(TC.gen.l):dim(TC.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  TC.Q.ext.pr[k] <- length(which(TC.Q.ext.sum > 0)) / iter
  print("THYLACOLEO")
  
  
  
  ## THYLACINUS (TH)
  ## set storage matrices & vectors
  TH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    TH.popmat <- TH.popmat.orig
    
    TH.n.mat <- matrix(0, nrow=TH.age.max+1,ncol=(t+1))
    TH.n.mat[,1] <- TH.init.vec

    for (i in 1:t) {
      # stochastic survival values
      TH.s.alpha <- estBetaParams(TH.Sx, TH.s.sd.vec^2)$alpha
      TH.s.beta <- estBetaParams(TH.Sx, TH.s.sd.vec^2)$beta
      TH.s.stoch <- rbeta(length(TH.s.alpha), TH.s.alpha, TH.s.beta)
      TH.s.stoch[1:TH.alpha] <- TH.s.stoch[1:TH.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/TH.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        TH.s.stoch <- TH.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      TH.fert.stch <- rnorm(length(TH.popmat[,1]), TH.pred.p.mm, TH.m.sd.vec)
      
      TH.totN.i <- sum(TH.n.mat[,i], na.rm=T)
      TH.pred.red <- TH.a.lp/(1+(TH.totN.i/TH.b.lp)^TH.c.lp)
      
      diag(TH.popmat[2:(TH.age.max+1),]) <- (TH.s.stoch[-(TH.age.max+1)])*TH.pred.red
      TH.popmat[TH.age.max+1,TH.age.max+1] <- 0 # (TH.s.stoch[TH.age.max+1])*TH.pred.red
      TH.popmat[1,] <- ifelse(TH.fert.stch < 0, 0, TH.fert.stch)
      TH.n.mat[,i+1] <- TH.popmat %*% TH.n.mat[,i]

    } # end i loop
    
    TH.n.sums.mat[e,] <- (as.vector(colSums(TH.n.mat)))

  } # end e loop
  
  # total N
  TH.n.md <- apply(TH.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  TH.n.up <- apply(TH.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  TH.n.lo <- apply(TH.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  TH.Q.ext.mat <- ifelse(TH.n.sums.mat < Q.ext.thresh, 1, 0)
  TH.Q.ext.sum <- apply(TH.Q.ext.mat[,ceiling(TH.gen.l):dim(TH.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  TH.Q.ext.pr[k] <- length(which(TH.Q.ext.sum > 0)) / iter
  print("THYLACINUS")
  

    
  ## SARCOPHILUS (SH)
  ## set storage matrices & vectors
  SH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    SH.popmat <- SH.popmat.orig
    
    SH.n.mat <- matrix(0, nrow=SH.age.max+1,ncol=(t+1))
    SH.n.mat[,1] <- SH.init.vec

    for (i in 1:t) {
      # stochastic survival values
      SH.s.alpha <- estBetaParams(SH.Sx, SH.s.sd.vec^2)$alpha
      SH.s.beta <- estBetaParams(SH.Sx, SH.s.sd.vec^2)$beta
      SH.s.stoch <- rbeta(length(SH.s.alpha), SH.s.alpha, SH.s.beta)
      SH.s.stoch[1:SH.alpha] <- SH.s.stoch[1:SH.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/SH.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        SH.s.stoch <- SH.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      SH.fert.stch <- rnorm(length(SH.popmat[,1]), SH.m.vec, SH.m.sd.vec)
      
      SH.totN.i <- sum(SH.n.mat[,i], na.rm=T)
      SH.pred.red <- SH.a.lp/(1+(SH.totN.i/SH.b.lp)^SH.c.lp)
  
      diag(SH.popmat[2:(SH.age.max+1),]) <- (SH.s.stoch[-1])*SH.pred.red
      SH.popmat[SH.age.max+1,SH.age.max+1] <- 0 # (SH.s.stoch[SH.age.max+1])*SH.pred.red
      SH.popmat[1,] <- ifelse(SH.fert.stch < 0, 0, SH.fert.stch)
      SH.n.mat[,i+1] <- SH.popmat %*% SH.n.mat[,i]

    } # end i loop
    
    SH.n.sums.mat[e,] <- (as.vector(colSums(SH.n.mat)))

  } # end e loop
  
  # total N
  SH.n.md <- apply(SH.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  SH.n.up <- apply(SH.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  SH.n.lo <- apply(SH.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  SH.Q.ext.mat <- ifelse(SH.n.sums.mat < Q.ext.thresh, 1, 0)
  SH.Q.ext.sum <- apply(SH.Q.ext.mat[,ceiling(SH.gen.l):dim(SH.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  SH.Q.ext.pr[k] <- length(which(SH.Q.ext.sum > 0)) / iter
  print("SARCOPHILUS")
  
  
  
  ## DASYURUS (DM)
  ## set storage matrices & vectors
  DM.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    DM.popmat <- DM.popmat.orig
    
    DM.n.mat <- matrix(0, nrow=DM.age.max+1,ncol=(t+1))
    DM.n.mat[,1] <- DM.init.vec

    for (i in 1:t) {
      # stochastic survival values
      DM.s.alpha <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$alpha
      DM.s.beta <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$beta
      DM.s.stoch <- rbeta(length(DM.s.alpha), DM.s.alpha, DM.s.beta)
      DM.s.stoch[1:DM.alpha] <- DM.s.stoch[1:DM.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/DM.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        DM.s.stoch <- DM.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      DM.fert.stch <- rnorm(length(DM.popmat[,1]), DM.m.vec, DM.m.sd.vec)
      
      DM.totN.i <- sum(DM.n.mat[,i], na.rm=T)
      DM.pred.red <- DM.a.lp/(1+(DM.totN.i/DM.b.lp)^DM.c.lp)
      
      diag(DM.popmat[2:(DM.age.max+1),]) <- (DM.s.stoch[-1])*DM.pred.red
      DM.popmat[DM.age.max+1,DM.age.max+1] <- 0 # (DM.s.stoch[DM.age.max+1])*DM.pred.red
      DM.popmat[1,] <- ifelse(DM.fert.stch < 0, 0, DM.fert.stch)
      DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]

    } # end i loop
    
    DM.n.sums.mat[e,] <- (as.vector(colSums(DM.n.mat)))

  } # end e loop
  
  # total N
  DM.n.md <- apply(DM.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  DM.n.up <- apply(DM.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  DM.n.lo <- apply(DM.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  DM.Q.ext.mat <- ifelse(DM.n.sums.mat < Q.ext.thresh, 1, 0)
  DM.Q.ext.sum <- apply(DM.Q.ext.mat[,ceiling(DM.gen.l):dim(DM.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  DM.Q.ext.pr[k] <- length(which(DM.Q.ext.sum > 0)) / iter
  print("DASYURUS")
  
  
  
  ## TACHYGLOSSUS (TA)
  ## set storage matrices & vectors
  TA.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    TA.popmat <- TA.popmat.orig
    
    TA.n.mat <- matrix(0, nrow=TA.age.max+1,ncol=(t+1))
    TA.n.mat[,1] <- TA.init.vec

    for (i in 1:t) {
      # stochastic survival values
      TA.s.alpha <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$alpha
      TA.s.beta <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$beta
      TA.s.stoch <- rbeta(length(TA.s.alpha), TA.s.alpha, TA.s.beta)
      TA.s.stoch[1:TA.alpha] <- TA.s.stoch[1:TA.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/TA.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        TA.s.stoch <- TA.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      TA.fert.stch <- rnorm(length(TA.popmat[,1]), TA.pred.p.mm, TA.m.sd.vec)
      
      TA.totN.i <- sum(TA.n.mat[,i], na.rm=T)
      TA.pred.red <- TA.a.lp/(1+(TA.totN.i/TA.b.lp)^TA.c.lp)
      
      diag(TA.popmat[2:(TA.age.max+1),]) <- (TA.s.stoch[-(TA.age.max+1)])*TA.pred.red
      TA.popmat[TA.age.max+1,TA.age.max+1] <- (TA.s.stoch[TA.age.max+1])*TA.pred.red
      TA.popmat[1,] <- ifelse(TA.fert.stch < 0, 0, TA.fert.stch)
      TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
      
    } # end i loop
    
    TA.n.sums.mat[e,] <- (as.vector(colSums(TA.n.mat)))

  } # end e loop
  
  # total N
  TA.n.md <- apply(TA.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  TA.n.up <- apply(TA.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  TA.n.lo <- apply(TA.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

  # quasi-extinction probability
  TA.Q.ext.mat <- ifelse(TA.n.sums.mat < Q.ext.thresh, 1, 0)
  TA.Q.ext.sum <- apply(TA.Q.ext.mat[,ceiling(TA.gen.l):dim(TA.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  TA.Q.ext.pr[k] <- length(which(TA.Q.ext.sum > 0)) / iter
  print("TACHYGLOSSUS")
  
  
  
  ## MEGALIBGWILIA (MR)
  ## set storage matrices & vectors
  MR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    MR.popmat <- MR.popmat.orig
    
    MR.n.mat <- matrix(0, nrow=MR.age.max+1,ncol=(t+1))
    MR.n.mat[,1] <- MR.init.vec

    for (i in 1:t) {
      # stochastic survival values
      MR.s.alpha <- estBetaParams(MR.Sx, MR.s.sd.vec^2)$alpha
      MR.s.beta <- estBetaParams(MR.Sx, MR.s.sd.vec^2)$beta
      MR.s.stoch <- rbeta(length(MR.s.alpha), MR.s.alpha, MR.s.beta)
      MR.s.stoch[1:MR.alpha] <- MR.s.stoch[1:MR.alpha] * (1 - juv.mort.inc.vec[k])
      
      if (rbinom(1, 1, 0.14/MR.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        MR.s.stoch <- MR.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      MR.fert.stch <- rnorm(length(MR.popmat[,1]), MR.pred.p.mm, MR.m.sd.vec)
      
      MR.totN.i <- sum(MR.n.mat[,i], na.rm=T)
      MR.pred.red <- MR.a.lp/(1+(MR.totN.i/MR.b.lp)^MR.c.lp)
      
      diag(MR.popmat[2:(MR.age.max+1),]) <- (MR.s.stoch[-(MR.age.max+1)])*MR.pred.red
      MR.popmat[MR.age.max+1,MR.age.max+1] <- (MR.s.stoch[MR.age.max+1])*MR.pred.red
      MR.popmat[1,] <- ifelse(MR.fert.stch < 0, 0, MR.fert.stch)
      MR.n.mat[,i+1] <- MR.popmat %*% MR.n.mat[,i]

    } # end i loop
    
    MR.n.sums.mat[e,] <- (as.vector(colSums(MR.n.mat)))

  } # end e loop
  
  # total N
  MR.n.md <- apply(MR.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  MR.n.up <- apply(MR.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  MR.n.lo <- apply(MR.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  MR.Q.ext.mat <- ifelse(MR.n.sums.mat < Q.ext.thresh, 1, 0)
  MR.Q.ext.sum <- apply(MR.Q.ext.mat[,ceiling(MR.gen.l):dim(MR.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  MR.Q.ext.pr[k] <- length(which(MR.Q.ext.sum > 0)) / iter
  print("MEGALIBGWILIA")
  
  
  print(juv.mort.inc.vec[k]) 
} # end k loop


## compile species outputs
spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)

DP.auc <- sum(sum(DP.Q.ext.pr)/length(juv.mort.inc.vec))
PA.auc <- sum(sum(PA.Q.ext.pr)/length(juv.mort.inc.vec))
ZT.auc <- sum(sum(ZT.Q.ext.pr)/length(juv.mort.inc.vec))
PH.auc <- sum(sum(PH.Q.ext.pr)/length(juv.mort.inc.vec))
VU.auc <- sum(sum(VU.Q.ext.pr)/length(juv.mort.inc.vec))
PG.auc <- sum(sum(PG.Q.ext.pr)/length(juv.mort.inc.vec))
SS.auc <- sum(sum(SS.Q.ext.pr)/length(juv.mort.inc.vec))
PT.auc <- sum(sum(PT.Q.ext.pr)/length(juv.mort.inc.vec))
SO.auc <- sum(sum(SO.Q.ext.pr)/length(juv.mort.inc.vec))
MN.auc <- sum(sum(MN.Q.ext.pr)/length(juv.mort.inc.vec))
OR.auc <- sum(sum(OR.Q.ext.pr)/length(juv.mort.inc.vec))
NR.auc <- sum(sum(NR.Q.ext.pr)/length(juv.mort.inc.vec))
GN.auc <- sum(sum(GN.Q.ext.pr)/length(juv.mort.inc.vec))
DN.auc <- sum(sum(DN.Q.ext.pr)/length(juv.mort.inc.vec))
AL.auc <- sum(sum(AL.Q.ext.pr)/length(juv.mort.inc.vec))
TC.auc <- sum(sum(TC.Q.ext.pr)/length(juv.mort.inc.vec))
TH.auc <- sum(sum(TH.Q.ext.pr)/length(juv.mort.inc.vec))
SH.auc <- sum(sum(SH.Q.ext.pr)/length(juv.mort.inc.vec))
DM.auc <- sum(sum(DM.Q.ext.pr)/length(juv.mort.inc.vec))
TA.auc <- sum(sum(TA.Q.ext.pr)/length(juv.mort.inc.vec))
MR.auc <- sum(sum(MR.Q.ext.pr)/length(juv.mort.inc.vec))

Qext.auc <- c(DP.auc,PA.auc,ZT.auc,PH.auc,VU.auc,PG.auc,SS.auc,PT.auc,SO.auc,MN.auc,OR.auc,NR.auc,GN.auc,DN.auc,AL.auc,TC.auc,TH.auc,SH.auc,DM.auc,TA.auc,MR.auc)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")
Qextpr.dat <- data.frame(juv.mort.inc.vec,DP.Q.ext.pr,PA.Q.ext.pr,ZT.Q.ext.pr,PH.Q.ext.pr,VU.Q.ext.pr,PG.Q.ext.pr,SS.Q.ext.pr,PT.Q.ext.pr,SO.Q.ext.pr,MN.Q.ext.pr,OR.Q.ext.pr,NR.Q.ext.pr,GN.Q.ext.pr,DN.Q.ext.pr,AL.Q.ext.pr,TC.Q.ext.pr,TH.Q.ext.pr,SH.Q.ext.pr,DM.Q.ext.pr,TA.Q.ext.pr,MR.Q.ext.pr)
colnames(Qextpr.dat) <- c("JSred",labs.vec)
dat.out <- data.frame(spp.mass.vec, Qext.auc)
colnames(dat.out) <- c("M","AUC")
rownames(dat.out) <- labs.vec

write.table(Qextpr.dat, file="juvsurvredQpr.csv", sep=",", dec = ".", row.names = F, col.names = T)
write.table(dat.out, file="juvsurvredout.csv", sep=",", dec = ".", row.names = T, col.names = T)

save.image("juvsurvred.RData")

