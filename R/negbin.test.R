#!/bin/env
## R --vanilla << "EOF" # Pipe all subsequent lines into R.

library(MASS)
library(lattice)
library(multcomp)
library(rjson)
library(DESeq)
#library(limma)

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

design_file = "design.txt"
contrast = "contrast.txt"
data_file = "data.txt"

data = read.table(data_file, header=T, row.names=1, sep=",")
features = rownames(data); nfeat = length(features)
samples = colnames(data); nsamples = length(samples)
groups = unique(unlist(lapply(strsplit(samples,".",fixed=T), "[[", 1))); ngroups = length(groups)

## Design matrix ##
design = read.table(design_file, header=T, row.names=1,  sep=",")
design = t(design)
covar = colnames(design)
ncovar = length(covar)
lvls = c()
for (cov in covar){ lvls = c(lvls,paste(cov,levels(as.factor(design[,cov])),sep=""))}
nlvls = length(lvls)

X = matrix(0,nsamples,length(lvls))
colnames(X) = lvls
rownames(X) = samples
for (i in 1:nsamples){
    gr = strsplit(samples[i],".",fixed=T)[[1]][1]
    for (j in 1:ncovar){
        cov = covar[j]
        val = design[gr,cov]
        X[samples[i],paste(cov,val,sep="")] = 1
    }
}

## Build the right part of the regression formula ##
formule = lvls[1]
for (c in lvls[2:nlvls]){ formule = paste(formule,"+",as.name(c)) }

## Calculate the model for each feature ##
results = list()
#for (i in 1:nrow(data)){
    i=1
    f = features[i]
    Y = t(data[f,])
    g = as.data.frame(cbind(Y,X))
    regressionFormula = formula(paste(features[i],"~",formule))

    nbmodel = glm.nb(regressionFormula, data=g) # AIC must be minimal
    summ = summary.glm(nbmodel)
    coeff = as.data.frame(summ$coefficients)

    result = matrix(NA,nlvls+1,4)
    rownames(result) = c("(Intercept)",lvls)
    colnames(result) = c("Estimate","Std. Error","t value","Pr(>|t|)")
    for (lvl in rownames(coeff)){
        for (res in colnames(coeff)){
            result[lvl,res] = coeff[lvl,res]
        }
    }
    results[[f]] = result

  ## Contrasts ##
  contrast = contrMat(rep(nfeat,ngroups), type="Tukey") # or Dunnett

  #test = glht(nbmodel, linfct=mcp(temp=contrast))
#}


