#!/bin/env
## R --vanilla << "EOF" # Pipe all subsequent lines into R.

library(MASS)
library(lattice)
library(multcomp)
library(rjson)

args=commandArgs(trailingOnly = TRUE)
args
set.seed(123)

## Set a fake experiment dataset
n = 30
features = paste(rep("feat",n), seq(n), sep="")
g1.1 = rnegbin(n, 200,2)
g1.2 = rnegbin(n, 200,2)
g1.3 = rnegbin(n, 200,2)
g2.1 = rnegbin(n, 100,0.2)
g2.2 = rnegbin(n, 100,0.2)
g2.3 = rnegbin(n, 100,0.2)
g3.1 = rnegbin(n, 700,0.4)
g3.2 = rnegbin(n, 700,0.4)
g3.3 = rnegbin(n, 700,0.4)
data = data.frame(g1.1,g1.2,g1.3,g2.1,g2.2,g2.3,g3.1,g3.2,g3.3, row.names=features)
write.table(data,"data.txt", sep=",", row.names=T, col.names=T, quote=F)


## Design matrix
design = "design.txt"
contrast = "contrast.txt"
filename = "data.txt"

#main <- function(data, design, contrast){

  data <- read.table(filename, header=T, row.names=1, sep=",")
  design <- read.table(design, header=T, row.names=1,  sep=",")
  design = as.data.frame(t(design))

  features = rownames(data)
  nfeat = length(features)
  cond = rownames(design)
  cov = colnames(design)
  ncond = length(cond)
  ncov = length(cov)

  # Build the right part of the regression formula
  f = cov[1]
  for (c in cov[2:ncov]){ f = paste(f,"+",c) }

  for (i in 1:length(design)){ design[,i] = as.factor(design[,i]) }
  estimate = matrix(,nfeat,ncov+1)
  stderror = matrix(,nfeat,ncov+1)
  zvalue = matrix(,nfeat,ncov+1)
  pvalue = matrix(,nfeat,ncov+1)

  for (i in 1:length(data[,1])){ # for each feature
    Y = as.data.frame(t(data[i,]))
    g = cbind(Y,design)
    ff = formula(paste(colnames(Y)[1],"~",f))

    nbmodel = glm.nb(ff, data=g)
    summ = summary(nbmodel)
    coef = as.data.frame(summ$coefficients)
    estimate[i,] = coef$"Estimate"
    stderror[i,] = coef$"Std. Error"
    zvalue[i,] = coef$"z value"
    pvalue[i,] = coef$"Pr(>|z|)"
  }

  data$estimate = estimate
  data$stderror = stderror
  data$zvalue = zvalue
  data$pvalue = pvalue

  ##contrast = glht(nbmodel, linfct=mcp(temp="Tukey"))

  ## contrast <- read.table(contrast.file, header=F, sep="\t")
  ## contrast.matrix <- data.matrix(contrast)

  #result = data.frame()
  #write.table(result, "ouput")
#  }


