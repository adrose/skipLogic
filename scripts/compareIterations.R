library(psych)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)
nsim1 <- rep(1000,8)      # how many simulations?
ne1 <- rep(1000,8)        # how many examinees per simulation?
nitems1 <- rep(20,8)      # how many items (total)?
probes1 <- c(rep(2,4), rep(6,4))       # how many items are probes?
amin1 <- rep(0.5,8)       # minimum value for a
amax1 <- rep(3.0,8)       # maximum value for a
dminp1 <- rep(-3,8)       # minimum value for d (probes)
dmaxp1 <- c(-2,-1,0,1,-2,-1,0,1)       # maximum value for d (probes)
dmin1 <- rep(0,8)         # minimum value for d
dmax1 <- rep(3,8)         # maximum value for d
cor_crit1 <- rep(0.8,8)


for(i in 1:8){
    nsim <- nsim1[i]      # how many simulations?
    ne <- ne1[i]        # how many examinees per simulation?
    nitems <- 20      # how many items (total)?
    probes <- probes1[i]       # how many items are probes?
    amin <- 0.5       # minimum value for a
    amax <- 3.0       # maximum value for a
    dminp <- -3       # minimum value for d (probes)
    dmaxp <- dmaxp1[i]       # maximum value for d (probes)
    dmin <- 0         # minimum value for d
    dmax <- 3         # maximum value for d
    cor_crit <- 0.8   # true correlation between criterion and theta
    
    v <- sqrt((1/cor_crit^2)-1)  # convert cor_crit to variance of residuals
    cl <- makeCluster(6)
    registerDoParallel(cl)
    results <- foreach(z=seq(1,nsim), .combine=rbind) %dopar% {
        library(psych)
        a <- runif(nitems,amin,amax)                                         # set item discriminations
        d <- c(runif(probes,dminp,dmaxp),runif((nitems-probes),dmin,dmax))   # set item difficulties
        dat <- sim.irt(nvar = nitems, n = ne, a=a, d=d)                      # simulate data
        try(modp <- irt.fa(dat$items,plot=FALSE))                            # calibrate item parameters before applying skip logic
        try(IRT_scores_no_SL <- scoreIrt(modp,dat$items)$theta1)             # calculate IRT scores before applying skip logic
        sum_no_SL <- rowSums(dat$items)                                      # calculate sum scores before applying skip logic
        for (i in 1:dim(dat$items)[1]) {if (rowSums(dat$items[,1:probes])[i] == 0) {dat$items[i,(probes+1):nitems] <- 0}}   # apply the skip logic
        try(mod <- irt.fa(dat$items,plot=FALSE))                             # calibrate item parameters after applying skip logic
        DV <- dat$theta + rnorm(ne,0,v)                                      # create criterion variable (DV to be predicted by scores)
        try(IRT_scores <- scoreIrt(mod,dat$items)$theta1)                    # calculate IRT scores after applying skip logic
        sum_scores <- rowSums(dat$items)                                     # calculate sum scores after applying skip logic
        x <- data.frame(DV,IRT_scores_no_SL,sum_no_SL,IRT_scores,sum_scores)
        cor(x)[2:5,1]
    }
    # Kill the cluster
    stopCluster(cl)
    
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
    
    write.csv(results,paste("sims",nsim,ne,nitems,probes,amin,amax,dminp,dmaxp,dmin,dmax,"final.csv", sep='_'))
}

# Now create a plot with the average corellation values from each of these output csv's
outCorVals <- NULL
for(i in 1:8){
    nsim <- nsim1[i]      # how many simulations?
    ne <- ne1[i]        # how many examinees per simulation?
    nitems <- 20      # how many items (total)?
    probes <- probes1[i]       # how many items are probes?
    amin <- 0.5       # minimum value for a
    amax <- 3.0       # maximum value for a
    dminp <- -3       # minimum value for d (probes)
    dmaxp <- dmaxp1[i]       # maximum value for d (probes)
    dmin <- 0         # minimum value for d
    dmax <- 3         # maximum value for d
    cor_crit <- 0.8   # true correlation between criterion and theta
    inputFile <- read.csv(paste("sims",nsim,ne,nitems,probes,amin,amax,dminp,dmaxp,dmin,dmax,"final.csv", sep='_'))
    outRow <- c(probes, dmaxp, mean(inputFile$IRT_no_SL),  mean(inputFile$IRT_after_SL))
    outCorVals <- rbind(outCorVals, outRow)
}
rownames(outCorVals) <- NULL
outCorVals <- as.data.frame(outCorVals)
## Now plot these values
outPlot <- ggplot(outCorVals, aes(x=V2)) +
  geom_point(aes(y=V3),color='blue', size=5) +
  geom_point(aes(y=V4),color='red', size=5) +
  facet_grid(V1 ~.) +
  xlab("dmaxp1") +
  ylab("Mean Cor Value")
pdf('meanIRTValueWandWOSL.pdf')
outPlot
dev.off()
