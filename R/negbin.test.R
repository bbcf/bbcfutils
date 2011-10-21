#!/bin/env
## R --vanilla << "EOF" # Pipe all subsequent lines into R.

library(MASS)
library(lattice)
library(multcomp)
library(limma)
library(rjson)

args=commandArgs(trailingOnly = TRUE)
args
set.seed(123)

## Create a fake experiment dataset
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

#main <- function(filename, design, contrast){

  data <- read.table(filename, header=T, row.names=1, sep=",")
  design <- read.table(design, header=T, row.names=1,  sep=",")
  design = as.data.frame(t(design))

  features = rownames(data)
  samples = rownames(design)
  covar = colnames(design) # covariates

  nfeat = length(features)
  nsamples = length(samples)
  ncovar = length(covar)

  # Determine groups
  groups = unique(unlist(lapply(strsplit(samples,".",fixed=T), "[[", 1)))

  # Build the right part of the regression formula
  f = covar[1]
  for (c in covar[2:ncovar]){ f = paste(f,"+",as.name(c)) }

  for (i in 1:ncovar){ design[,i] = as.factor(design[,i]) }
  estimate = matrix(,nfeat,ncovar+1)
  stderror = matrix(,nfeat,ncovar+1)
  zvalue = matrix(,nfeat,ncovar+1)
  pvalue = matrix(,nfeat,ncovar+1)

  for (i in 1:1){#nrow(data)){ # for each feature
    Y = as.data.frame(t(data[i,]))
    g = cbind(Y,design)
    regressionFormula = formula(paste(colnames(Y)[1],"~",f))

    nbmodel = glm.nb(regressionFormula, data=g)
    summ = summary(nbmodel)
    coeff = as.data.frame(summ$coefficients)
    estimate[i,] = coeff$"Estimate" # 'beta' coefficients of the regression
    stderror[i,] = coeff$"Std. Error"
    zvalue[i,] = coeff$"z value"
    pvalue[i,] = coeff$"Pr(>|z|)"

    contrast = contrMat(rep(nfeat,nsamples), type="Tukey")
    test = glht(nbmodel, linfct=mcp(temp=contrast))
    #contrast = glht(nbmodel, linfct=mcp(temp="Tukey"))
    ## contrast <- read.table(contrast.file, header=F, sep="\t")
    ## contrast.matrix <- data.matrix(contrast)

  }

  colnames(estimate) <- paste(rep("estimate",3),".",groups,sep="")
  colnames(stderror) <- paste(rep("stderror",3),".",groups,sep="")
  colnames(pvalue) <- paste(rep("pvalue",3),".",groups,sep="")
  estimate = as.data.frame(estimate)
  stderror = as.data.frame(stderror)
  pvalue = as.data.frame(pvalue)
  data = cbind(data, estimate, stderror, pvalue)

  #result = data.frame()
  #write.table(result, "ouput")
#  }


