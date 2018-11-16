library(psych)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)
library(caret)
library(mirt)

## Load the data
ocd.data <- read.csv('./OCD.csv')

## Now run all of the items through MIRT
# First remove all completley empty rows
ocd.data <- ocd.data[-which(apply(ocd.data, 1, function(x) sum(is.na(x)))==17),]
mod1 <- mirt(data=ocd.data[,c(2:18)], 1, IRTpars=T)


## Now obtain our difficulty estimates
outVals <- coef(mod1, irtPars=T)
