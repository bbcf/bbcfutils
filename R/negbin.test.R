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
means1 = sample(100:200,n,replace=T); thetas1 = sample(2:5,n,replace=T)/10
means2 = sample(100:200,n,replace=T); thetas2 = sample(2:5,n,replace=T)/10
means3 = sample(500:700,n,replace=T); thetas3 = sample(8:12,n,replace=T)/10
features = paste(rep("feat",n), seq(n), sep="")
data = data.frame(row.names=samples)
for (i in 1:n){
    data[features[i]] = c(rnegbin(3,means1[i],thetas1[i]),rnegbin(3,means2[i],thetas2[i]),rnegbin(3,means3[i],thetas3[i]))
}
data = t(data)
write.table(data,"tests/data.txt", sep=",", row.names=T, col.names=T, quote=F)


main <- function(filename, design, contrast){

    design_file = "tests/design.txt"
    contrast_file = "tests/contrast.txt"
    data_file = "tests/data.txt"

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
    for (cov in covar){
        lvls = c(lvls,paste(cov,levels(as.factor(design[,cov])),sep=""))
    }
    nlvls = length(lvls)

    ## Contrasts matrix ##
    if (file.exists(contrast_file)){
        contrast = as.matrix(read.table(contrast_file, header=T, row.names=1, sep=","))
    }else{
        contrast = contrMat(rep(nfeat,ngroups), type="Tukey") # or Dunnett
    }
    ncomp = length(rownames(contrast))

    ## Regression coefficients
    X = c()
    for (j in 1:nsamples){
      gr = strsplit(samples[j],".",fixed=T)[[1]][1]
      i = which(groups==gr)
      e = paste(covar[1],design[i,1], sep="")
      for (j in 2:ncovar){
         e = paste(e,paste(covar[j],design[i,j], sep=""), sep="_")
      }
      X = cbind(X,e)
    }

    ## Calculate the model for each feature ##
    results = list()
    comparisons = list()
    for (i in 1:nrow(data)){
        f = features[i]
        Y = t(data[f,])
        g = data.frame(Y=Y,All=as.factor(t(X)))
        regressionFormula = formula(paste(features[i],"~",as.name("All")))

        nbmodel = glm.nb(regressionFormula, data=g) # AIC must be minimal
        nbmodel.summ = summary.glm(nbmodel)
        coeff = as.data.frame(nbmodel.summ$coefficients)
        lvls = rownames(coeff)
        for (i in 1:length(lvls)){
            if (lvls[i] != "(Intercept)"){ lvls[i] = strsplit(lvls[i],"All")[[1]][2]}
        }
        rownames(coeff) = lvls

        result = matrix(NA,ncovar,4)
        rownames(result) = rownames(coeff)
        colnames(result) = colnames(coeff)
        for (lvl in lvls){
            for (res in colnames(coeff)){
                result[lvl,res] = coeff[lvl,res]
            }
        }
        results[[f]] = result

        ## Contrasts ##
        comp = glht(nbmodel, linfct=mcp(All=contrast))
        comp.summ = summary(comp)
        estimate = comp.summ$test$coefficients
        pval = comp.summ$test$pvalues

        comparisons[[f]] = data.frame("Estimate"=estimate, "Pvalue"=pval)
        warns = warnings()
    }

    ## Compute adjusted p-values ##
    pvalues = c()
    for (i in 1:length(rownames(contrast))){
        pvalues = c(pvalues,unlist(lapply(lapply(comparisons,"[[",2),"[[",i)))
    }
    pvalues.adj = p.adjust(pvalues,"BH")
    for (i in 1:nfeat){
        comparisons[[features[i]]]["Adj.Pvalue"] = pvalues.adj[nfeat*(0:(ncomp-1))+i]
    }

    ## Return ##
    comparisons
}

comparisons = main(data_file, design_file, contrast_file)
unlink("negbin.test.txt") #deletes the file if it already exists
for (f in features) {
    write.table(comparisons[f],"negbin.test.txt",quote=F,row.names=T,col.names=T,append=T,sep="\t")
    write(c(),"negbin.test.txt",append=T)
}
