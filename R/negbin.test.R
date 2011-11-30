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


main <- function(data_file, design_file, contrast_file, nsamples, sep){

    data = read.table(data_file, header=T, row.names=1, sep=sep)
    data = data[,1:nsamples]
    features = rownames(data); nfeat = length(features)
    samples = colnames(data)
    groups = unique(unlist(lapply(strsplit(samples,".",fixed=T), "[[", 2))); ngroups = length(groups)

    ## Design matrix ##
    design = read.table(design_file, header=T, row.names=1,  sep=sep)
    design = t(design)
    covar = colnames(design)
    ncovar = length(covar)

    ## Contrasts matrix ##
    if (file.exists(contrast_file)){
        contrast = as.matrix(read.table(contrast_file, header=T, row.names=1, sep=sep))
    }else{
        contrast = contrMat(rep(nfeat,ngroups), type="Tukey") # or Dunnett
    }
    contrast.names = rownames(contrast)
    ncomp = length(contrast.names)

    ## Regression coefficients
    X = c()
    for (j in 1:nsamples){
      gr = strsplit(samples[j],".",fixed=T)[[1]][2]
      i = which(groups==gr)
      e = paste(covar[1],design[i,1], sep="")
      for (j in 2:ncovar){
         e = paste(e,paste(covar[j],design[i,j], sep=""), sep="_")
      }
      X = cbind(X,e)
    }

    ## Calculate the model for each feature ##
    models = list() #contains all nbmodels
    comparisons = list() #contains all contrast summaries
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

        result = matrix(NA,length(lvls),4)
        rownames(result) = rownames(coeff)
        colnames(result) = colnames(coeff)
        for (lvl in lvls){
            for (res in colnames(coeff)){
                result[lvl,res] = coeff[lvl,res]
            }
        }
        models[[f]] = result

        ## Contrasts ##
        comp = glht(nbmodel, linfct=mcp(All=contrast))
        comp.summ = summary(comp)
        estimate = comp.summ$test$coefficients
        pval = comp.summ$test$pvalues
        comparisons[[f]] = data.frame("Estimate"=estimate, "Pvalue"=pval)
    }

    ## Compute adjusted p-values ##
    pvalues = c()
    for (i in 1:length(contrast.names)){
        pvalues = c(pvalues,unlist(lapply(lapply(comparisons,"[[",2),"[[",i)))
    }
    pvalues.adj = p.adjust(pvalues,"BH")
    for (i in 1:nfeat){
        comparisons[[features[i]]]["Adj.Pvalue"] = pvalues.adj[nfeat*(0:(ncomp-1))+i]
    }

    ## Group by comparison type and sort ##
    bycomp = list()
    comparisons.t = lapply(comparisons,FUN=t)
    for (i in 1:ncomp){
        a = c()
        for (f in features){
            a = cbind(a,comparisons.t[[f]][,i])
        }
        b = as.data.frame(t(a))
        rownames(b) = features
        b = b[order(b[,3]),] # sort w.r.t adjusted p-values
        bycomp[[contrast.names[i]]] = b
    }

    bycomp
}


create_fake_dataset <- function(n){
    samples = c("g1.1","g1.2","g1.3","g2.1","g2.2","g2.3","g3.1","g3.2","g3.3")
    samples = c(paste("counts",samples,sep="."),paste("rpkm",samples,sep="."))
    means1 = sample(100:200,n,replace=T); thetas1 = sample(2:5,n,replace=T)/10
    means2 = sample(100:200,n,replace=T); thetas2 = sample(2:5,n,replace=T)/10
    means3 = sample(500:700,n,replace=T); thetas3 = sample(8:12,n,replace=T)/10
    features = paste(rep("feat",n), seq(n), sep="")
    data = data.frame(row.names=samples)
    for (i in 1:n){
        line = c(rnegbin(3,means1[i],thetas1[i]),rnegbin(3,means2[i],thetas2[i]),rnegbin(3,means3[i],thetas3[i]))
        data[features[i]] = c(line, line/3)
    }
    data = as.data.frame(signif(t(data),2))
    write.table(data,"tests/data.txt", sep=",", row.names=T, col.names=T, quote=F)
}


test0 <- function(){
    data_file = "tests/data.txt"
    design_file = "tests/design.txt"
    contrast_file = "tests/contrast.txt"
    comparisons = main(data_file, design_file, contrast_file, nsamples=9, sep=",")
    unlink("tests/negbin.test.txt") #deletes the file if it already exists
    for (c in names(comparisons)) {
        result = signif(comparisons[[c]],4)
        write(c,"tests/negbin.test.txt",append=T)
        write.table(result,"tests/negbin.test.txt",quote=F,row.names=T,col.names=T,append=T,sep="\t")
        write(c(),"tests/negbin.test.txt",append=T)
    }
}

test1 <- function(){
    #missing replicates
    data_file = "/archive/epfl/bbcf/jdelafon/MEF/mef_genes.csv"
    design_file = "tests/design_mef.txt"
    contrast_file = "tests/contrast_mef.txt"
    comparisons = main(data_file, design_file, contrast_file, nsamples=2, sep="\t")
}

test2 <- function(){
    data_file = "tests/mult_genes.csv"
    design_file = "tests/design_mef.txt"
    contrast_file = "tests/contrast_mef.txt"
    main(data_file, design_file, contrast_file, nsamples=6, sep="\t")
}
