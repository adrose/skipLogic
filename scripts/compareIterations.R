library(psych)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)
library(caret)
# First identify all of our conditions
probes <- c(2,4,6)
amin <- c(0.3,1.5)
dmin1 <- c(-1,1)
dmin <- c(-3,-1)
theta_dist <- c(1,2)
all.folds <- expand.grid(probes,amin,dmin,dmin1,theta_dist)

cl <- makeCluster(32)
registerDoParallel(cl)
tmp <- foreach(E=1:dim(all.folds)[1], .export='all.folds') %do% {
  library(psych)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(reshape2)
  library(caret)
  nsim <- 2000      # how many simulations?
  ne <- 10000        # how many examinees per simulation?
  nitems <- 20      # how many items (total)?
  probes <- all.folds[E,1]       # how many items are probes?
  amin <- all.folds[E,2]       # minimum value for a
  if (amin == 0.3) {amax <- 1.5}
  if (amin == 1.5) {amax <- 3.5}
  dminp <- all.folds[E,3]       # minimum value for d (probes)
  dmaxp <- dminp + 2       # maximum value for d (probes)
  dmin <- all.folds[E,4]         # minimum value for d
  dmax <- dmin + 2        # maximum value for d
  theta_dist <- all.folds[E,5]	# shape of theta function
  if (theta_dist == 1) {theta.sim <- rnorm(ne)}                            # set simulated theta
  if (theta_dist == 2) {theta.sim <- scale(rnorm(ne)^2)}
  cor_crit <- 0.8   # true correlation between criterion and theta   
  v <- sqrt((1/cor_crit^2)-1)  # convert cor_crit to variance of residuals
  # Now make a temporary output directory
  out.dir.name <- paste("sims",probes,amin,amax,dminp,dmaxp,dmin,dmax,theta_dist,E,sep=',')
  dir.create(out.dir.name, showWarnings = FALSE)
  results <- foreach(z=seq(1,nsim), .combine=rbind, .export='out.dir.name') %dopar% {
        if(!file.exists(paste(out.dir.name,"/", z, ".csv", sep=''))){
          library(psych)
          a <- runif(nitems,amin,amax)                                         # set item discriminations
          d <- c(runif(probes,dminp,dmaxp),runif((nitems-probes),dmin,dmax))   # set item difficulties
          dat <- sim.irt(nvar = nitems, n = ne, a=a, d=d, theta=theta.sim)                      # simulate data
          modp <- irt.fa(dat$items,plot=FALSE)                           # calibrate item parameters before applying skip logic
          IRT_scores_no_SL <- scoreIrt(modp,dat$items)$theta1             # calculate IRT scores before applying skip logic
          sum_no_SL <- rowSums(dat$items)                                      # calculate sum scores before applying skip logic
          for (i in 1:dim(dat$items)[1]) {if (rowSums(dat$items[,1:probes])[i] == 0) {dat$items[i,(probes+1):nitems] <- 0}}   # apply the skip logic
          mod <- irt.fa(dat$items,plot=FALSE)                             # calibrate item parameters after applying skip logic
          DV <- dat$theta + rnorm(ne,0,v)                                      # create criterion variable (DV to be predicted by scores)
          IRT_scores <- scoreIrt(mod,dat$items)$theta1                    # calculate IRT scores after applying skip logic
          sum_scores <- rowSums(dat$items)                                     # calculate sum scores after applying skip logic
          x <- data.frame(DV,IRT_scores_no_SL,sum_no_SL,IRT_scores,sum_scores)
          print(z)
          write.csv(cor(x)[2:5,1], paste(out.dir.name,"/", z, ".csv", sep=''), quote=F, row.names=F)
          cor(x)[2:5,1]
      }
      else{
        t(read.csv(paste(out.dir.name,"/", z, ".csv", sep='')))
      }
    }
    
    # Now create the DV variable
    a <- runif(nitems,amin,amax)                                         # set item discriminations
    d <- c(runif(probes,dminp,dmaxp),runif((nitems-probes),dmin,dmax))   # set item difficulties
    dat <- sim.irt(nvar = nitems, n = ne, a=a, d=d)                      # simulate data
    try(modp <- irt.fa(dat$items,plot=FALSE))                            # calibrate item parameters before applying skip logic
    try(IRT_scores_no_SL <- scoreIrt(modp,dat$items)$theta1)             # calculate IRT scores before applying skip logic
    sum_no_SL <- rowSums(dat$items)                                      # calculate sum scores before applying skip logic
    for (i in 1:dim(dat$items)[1]) {if (rowSums(dat$items[,1:probes])[i] == 0) {dat$items[i,(probes+1):nitems] <- 0}}   # apply the skip logic
    try(mod <- irt.fa(dat$items,plot=FALSE))                             # calibrate item parameters after applying skip logic
    DV <- dat$theta + rnorm(ne,0,v)
    
    results <- data.frame(results,nsim,ne,nitems,probes,amin,amax,dminp,dmaxp,dmin,dmax,cor(cbind(DV,dat$theta))[2,1])
    colnames(results)[1:4] <- c("IRT_no_SL","sum_score_no_SL","IRT_after_SL","sum_score_after_SL")
    
    write.csv(results,paste("sims",nsim,ne,nitems,probes,amin,amax,dminp,dmaxp,dmin,dmax,theta_dist,"final.csv", sep='_'))
}
# Kill the cluster
stopCluster(cl)
q()
