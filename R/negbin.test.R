#!/bin/env
# R --vanilla negbin.test.r --args -s '\t' -d design_file -c contrast_file -o output_file

library(MASS)
library(lattice)
library(multcomp)

args=commandArgs(trailingOnly = TRUE)
args
set.seed(123)

data_file = args[1]
sep = args[grep("-s",args)+1]
design_file = args[grep("-d",args)+1]
contrast_file = args[grep("-c",args)+1]
output_file = args[grep("-o",args)+1]


main <- function(data_file, sep="\t", contrast_file=FALSE, design_file=FALSE, output_file=FALSE){
    data = read.table(data_file, header=T, row.names=1, sep=sep)
    header = colnames(data)
    counts = grep("counts",header,fixed=T)
    nsamples = length(counts)
    data = data[,counts]

    ## Choose GLM if every group has replicates, DESeq otherwise ##
    a = c()
    for (g in groups){ if (length(which(conds==g))>1) {a = c(a,1)} }
    if (length(a) == length(groups) && file.exists(design_file) && file.exists(contrast_file)){
        print("GLM")
        design = read_design(design_file, sep)
        contrast = read_contrast(contrast_file, sep)
        comparisons = GLM(data, design, contrast, output_file)
    }else{
        print("DESeq")
        different = DES(data)
    }
}


read_design <- function(design_file, sep){
    design = read.table(design_file, header=T, row.names=1,  sep=sep)
    design = t(design)
}
read_contrast <- function(contrast_file, sep){
    contrast = as.matrix(read.table(contrast_file, header=T, row.names=1, sep=sep))
}


DES <- function(data){  ## DESeq ##
    library(DESeq)
    samples = colnames(data)
    conds = unlist(lapply(strsplit(samples,".",fixed=T), "[[", 2))

    result = list()
    cds <- newCountDataSet(data, conds)
    cds <- estimateSizeFactors(cds)
    cds <- estimateVarianceFunctions(cds)
    couples = combn(unique(groups),2)
    for (i in 1:dim(couples)[2]){
        res <- nbinomTest(cds, couples[1,i], couples[2,i])
        res = res[order(res[,8]),] # sort w.r.t. adjusted p-value
        result[[paste(couples[1,i],"-",couples[2,i])]] = res
    }
    ## Return ##
    result
}


GLM <- function(data, design, contrast, output_file=FALSE){
    features = rownames(data)
    samples = colnames(data)
    groups = unique(unlist(lapply(strsplit(samples,".",fixed=T), "[[", 2)))
    covar = colnames(design)
    contrast.names = rownames(contrast)

    nfeat = length(features); ngroups = length(groups); ncovar = length(covar)
    nsamples = length(samples); ncomp = length(contrast.names)

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
        b = b[order(b[,3]),] # sort w.r.t adjusted p-value
        bycomp[[contrast.names[i]]] = b
    }

    ## Return ##
    if (output_file == FALSE){
        print(bycomp)
    }else{
        for (comp in contrast.names) {
            result = signif(bycomp[[comp]],4)
            write(comp,output_file,append=T)
            write.table(result,output_file, quote=F,row.names=T,col.names=T,append=T,sep="\t")
            write(c(),output_file,append=T)
        }
    }
}


