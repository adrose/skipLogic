## This script will be used to plot all of the main effects
## of skip logic
## It'll be a beast of a figure
## FIrst load all library(s)
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

## Now go through and read all data, store all of the values in one mega csv
all.values <- foreach(E=1:dim(all.folds)[1], .combine='rbind') %do% {
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
  tmp <- read.csv(paste("sims",nsim,ne,nitems,probes,amin,amax,dminp,dmaxp,dmin,dmax,theta_dist,"final.csv", sep='_'))
  tmp$theta_dist <- theta_dist
  tmp
}
all.values$toPlot <- all.values$IRT_no_SL - all.values$IRT_after_SL

## Now go through all of oour main effects and find the differences
vals <- names(all.values)[c(9,12,14,17,10)]
plotVals <- NULL
for(i in vals){
  to.mod <- summarySE(data=all.values, measurevar='toPlot', groupvars=i)
  colnames(to.mod)[1] <- 'meanVals'
  to.mod$val <- i
  plotVals <- rbind(plotVals, to.mod)
}

## now create the plot
plotVals$val <- c(rep('Number of probe items', 3), rep("Probe item difficulties", 2), rep("Non-Probe item difficulties",2), rep("Thea Distribution", 2), rep("Item discrimination", 2))

## Now create the plot
outPlot <- ggplot(plotVals, aes(x=factor(meanVals), y=toPlot)) +
  geom_bar(stat='identity', position=position_dodge(), width=.5) + 
  geom_errorbar(aes(ymin=as.numeric(as.character(toPlot))-se, ymax=as.numeric(as.character(toPlot))+se), 
                       width = .5, position=position_dodge()) +
  theme_bw() +
  facet_grid(.~val, space = "free_x", scales='free_x') +
  ylab("Bias Due to Skip-Logic") +
  xlab("Conditions")

pdf("slSuppFig1.pdf")
outPlot
dev.off()
