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
samples = c("g1.1","g1.2","g1.3","g2.1","g2.2","g2.3","g3.1","g3.2","g3.3")
means1 = sample(100:200,30,replace=T); thetas1 = sample(1:10,30,replace=T)/10
means2 = sample(500:600,30,replace=T); thetas2 = sample(5:15,30,replace=T)/10
means3 = sample(700:900,30,replace=T); thetas3 = sample(10:20,30,replace=T)/10
features = paste(rep("feat",n), seq(n), sep="")
data = data.frame(row.names=samples)
for (i in 1:n){
    data[features[i]] = c(rnegbin(3,means1[i],thetas1[i]),rnegbin(3,means2[i],thetas2[i]),rnegbin(3,means3[i],thetas3[i]))
}
data = t(data)
write.table(data,"data.txt", sep=",", row.names=T, col.names=T, quote=F)

#main <- function(filename, design, contrast){

design = "design.txt"
contrast = "contrast.txt"
filename = "data.txt"

data <- read.table(filename, header=T, row.names=1, sep=",")
features = rownames(data); nfeat = length(features)
samples = colnames(data); nsamples = length(samples)
groups = unique(unlist(lapply(strsplit(samples,".",fixed=T), "[[", 1))); ngroups = length(groups)

## Design matrix ##
design = read.table(design, header=T, row.names=1,  sep=",")
design = as.data.frame(t(design))

## Covariates ##
covar = colnames(design)
ncovar = length(covar)

## Build the right part of the regression formula ##
formule = covar[1]
for (c in covar[2:ncovar]){ formule = paste(formule,"+",as.name(c)) }

## Initialization ##
for (i in 1:ncovar){ design[,i] = as.factor(design[,i]) }
estimate = matrix(,nfeat,ncovar)
stderror = matrix(,nfeat,ncovar)
zvalue = matrix(,nfeat,ncovar)
pvalue = matrix(,nfeat,ncovar)

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


