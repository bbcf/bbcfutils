#!/bin/env
## R --vanilla << "EOF" # Pipe all subsequent lines into R.

library(MASS)
library(lattice)
library(multcomp)
#library(limma)
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

#main <- function(filename, design, contrast){
  #contrast <- read.table(contrast, header=F, sep=",")
  #design <- read.table(design, header=T, row.names=1,  sep=",")
  #design = as.data.frame(t(design))

design = "design.txt"
contrast = "contrast.txt"
filename = "data.txt"

data <- read.table(filename, header=T, row.names=1, sep=",")
features = rownames(data); nfeat = length(features)
samples = colnames(data); nsamples = length(samples)
groups = unique(unlist(lapply(strsplit(samples,".",fixed=T), "[[", 1))); ngroups = length(groups)

## Design matrix ## (replicates as lines, groups as columns)
design = data.frame(row.names=groups)
for (s in 1:nsamples){
  newcol = rep(0,ngroups);
  for (g in 1:ngroups){
    if (strsplit(samples[s],'.',fixed=T)[[1]][1]==groups[g])
      {newcol[g]=1} 
    }
  design[samples[s]] = newcol
  }
design = t(design)
#design = cbind(data.frame(m=rep(1,nsamples)), design)

## Covariates ## (same as groups for us)
covar = colnames(design)
ncovar = length(covar)

## Build the right part of the regression formula ##
formule = covar[1]
for (c in covar[2:ncovar]){ formule = paste(formule,"+",as.name(c)) }

## Initialization ##
for (i in 1:ncovar){ design[,i] = as.factor(design[,i]) }
estimate = matrix(,nfeat,ncovar+1)
stderror = matrix(,nfeat,ncovar+1)
zvalue = matrix(,nfeat,ncovar+1)
pvalue = matrix(,nfeat,ncovar+1)

#for (i in 1:nrow(data)){
  i=1
  F = data[features[i],]
  Y = as.data.frame(t(F))
  g = cbind(Y,design)
  regressionFormula = formula(paste(features[i],"~",formule))

  nbmodel = glm.nb(regressionFormula, data=g) # AIC must be minimal
  summ = summary(nbmodel)
  coeff = as.data.frame(summ$coefficients)
  estimate[i,] = coeff$"Estimate" # 'beta' coefficients of the regression
  stderror[i,] = coeff$"Std. Error"
  zvalue[i,] = coeff$"z value"
  pvalue[i,] = coeff$"Pr(>|z|)"

  ## Contrasts ##
  contrast = contrMat(rep(nfeat,ngroups), type="Tukey") # ou Dunnett

  test = glht(nbmodel, linfct=mcp(temp=contrast))
  #contrast.matrix <- data.matrix(contrast)
#}

colnames(estimate) <- paste(rep("estimate",3),".",groups,sep="")
colnames(stderror) <- paste(rep("stderror",3),".",groups,sep="")
colnames(pvalue) <- paste(rep("pvalue",3),".",groups,sep="")
estimate = as.data.frame(estimate)
stderror = as.data.frame(stderror)
pvalue = as.data.frame(pvalue)
data = cbind(data, estimate, stderror, pvalue)

#result = data.frame()
#write.table(result, "ouput")


