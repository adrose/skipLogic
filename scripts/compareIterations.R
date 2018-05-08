library(psych)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)
library(caret)
# First identify all of our conditions
probes <- c(2,4,6)
amin <- c(0.3, 1.5)
dmin <- c(-2, -1, 0, 1)
dmin1 <- c(-2, -1, 0, 1)
all.folds <- expand.grid(probes,amin,dmin,dmin1)


for(i in 1:dim(all.folds)[1]){
    nsim <- 2000      # how many simulations?
    ne <- 1000        # how many examinees per simulation?
    nitems <- 20      # how many items (total)?
    probes <- all.folds[i,1]       # how many items are probes?
    amin <- all.folds[i,2]       # minimum value for a
    if(amin==0.3){amax <- 1.5}   # maximum value for a
    if(amin==1.5){amax <- 3.5}   # maximum value for a
    dminp <- all.folds[i,3]       # minimum value for d (probes)
    dmaxp <- dminp + 1       # maximum value for d (probes)
    dmin <- all.folds[i,4]         # minimum value for d
    dmax <- dmin + 1        # maximum value for d
    cor_crit <- 0.8   # true correlation between criterion and theta
    
    v <- sqrt((1/cor_crit^2)-1)  # convert cor_crit to variance of residuals
    cl <- makeCluster(8)
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
for(i in 1:dim(all.folds)[1]){
    nsim <- 2000      # how many simulations?
    ne <- 1000        # how many examinees per simulation?
    nitems <- 20      # how many items (total)?
    probes <- all.folds[i,1]       # how many items are probes?
    amin <- all.folds[i,2]       # minimum value for a
    if(amin==0.3){amax <- 1.5}   # maximum value for a
    if(amin==1.5){amax <- 3.5}   # maximum value for a
    dminp <- all.folds[i,3]       # minimum value for d (probes)
    dmaxp <- dminp + 1       # maximum value for d (probes)
    dmin <- all.folds[i,4]         # minimum value for d
    dmax <- dmin + 1        # maximum value for d
    cor_crit <- 0.8  # true correlation between criterion and theta
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
