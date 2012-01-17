#!usr/bin/env

# This script looks for differential expression in RNA-seq data
# of genomic features in different conditions.

# --------
#  Usage:
# --------
# R --slave -f negbin.test.R --args 'data_file' -s ',' -d 'design_file' -c 'contrast_file' -o 'output_file'

# Personal example:
# (On Vital-IT, .bashrc contains: alias R="bsub -qbbcf -Ip -XF R --vanilla")
# DESeq:
# R --slave -f negbin.test.R --args 'tests/testing_files/genes_expression.tab' -s 'tab' -o 'tests/testing_files/output_deseq.txt'
# GLM:
# R --slave -f negbin.test.R --args 'tests/testing_files/genes_expression.tab' -s 'tab' -d 'tests/testing_files/design_mado.txt' -c 'tests/testing_files/contrast_mado.txt' -o 'tests/testing_files/output_glm.txt'

# Arguments:
# -s sep: character separating fields in the input files. Use 'tab' for tab delimiter (not '\t').
# -d design: name of the file containing the design matrix (see below).
# -c contrast: name of the file containing the contrast matrix (see below).
# -o output_file: name of the file containing the results. It is either created, or appends the
#    results at the end of the file if it already exists.

# The script is made to take as input the files returned by rnaseq.py, which format is the following:
# - tab-delimited or CSV file containing (maybe normalized) read counts
# - lines representing genomic features
# - columns represent different samples
# - the first column contains feature names
# - columns containing counts are labeled "counts.<groupName>.<sampleIndex>" (e.g. counts.G1.3).

# If each group contains several replicate samples, a General Linear Model is fitted.
# If at least one group contains no replicates, DESeq is used.
#
# To fit the GLM, the user needs to provide a design matrix and a contrast matrix in text format
# (see tutorial [...]). If not found, DESeq is run.
# If no output_file is provided, output is printed to stdout.

# Still to implement:
# - Automatically generated contrast_file if not provided, just comparing groups (Tukey matrix -> ANOVA)


library(MASS)
library(multcomp)

args=commandArgs(trailingOnly = TRUE)
args

data_file = args[1]
sep = args[grep("-s",args)+1]
design_file = args[grep("-d",args)+1]
contrast_file = args[grep("-c",args)+1]
output_file = args[grep("-o",args)+1]
if (sep=='tab') sep='\t'


main <- function(data_file, sep="\t", contrast_file=FALSE, design_file=FALSE, output_file=FALSE){
    data = read.table(data_file, header=T, row.names=1, sep=sep)
    header = colnames(data)
    counts = grep("counts",header,fixed=T)
    data = round(data[,counts])

    ## Choose GLM if every group has replicates, DESeq otherwise ##
    samples = colnames(data)
    conds = unlist(lapply(strsplit(samples,".",fixed=T), "[[", 2))
    print(paste("- All groups have several runs: ", all(table(conds)>1)))
    print(paste("- Design file exists: ", file.exists(design_file)," - ",design_file))
    print(paste("- Contrast file exists: " ,file.exists(contrast_file)," - ",contrast_file))
    if (all(table(conds)>1) && file.exists(design_file) && file.exists(contrast_file)){
        print("Method: GLM")
        design = read_design(design_file, sep)
        contrast = read_contrast(contrast_file, sep)
        GLM(data, design, contrast, output_file)
    }else{
        print("Method: DESeq (either design/contrast files missing, or at least one group has no replicates.)")
        DES(data, output_file)
    }
}


read_design <- function(design_file, sep){
    design = read.table(design_file, header=T, row.names=1,  sep=sep)
    design = t(design)
}
read_contrast <- function(contrast_file, sep){
    contrast = as.matrix(read.table(contrast_file, header=T, row.names=1, sep=sep))
}


DES <- function(data, output_file=FALSE){  ## DESeq ##
    library(DESeq)
    samples = colnames(data)
    conds = unlist(lapply(strsplit(samples,".",fixed=T), "[[", 2))
    groups = unique(conds)

    result = list()
    cds <- newCountDataSet(data, conds)
    cds <- estimateSizeFactors(cds)
    cds <- estimateVarianceFunctions(cds, method='blind')
    couples = combn(unique(groups),2)
    contrast.names = c()
    for (i in 1:dim(couples)[2]){
        res <- nbinomTest(cds, couples[1,i], couples[2,i])
        res = res[order(res[,8]),] # sort w.r.t. adjusted p-value
        comp = paste(couples[1,i],"-",couples[2,i])
        contrast.names = c(contrast.names,comp)
        result[[comp]] = res
    }
    ## Return ##
    if (output_file == FALSE){
        print(result)
    }else{
        for (comp in contrast.names) {
            #result = signif(result[[comp]],4)
            write(comp,output_file,append=T)
            write.table(result[[comp]],output_file,quote=F,row.names=T,col.names=T,append=T,sep="\t")
            write(c(),output_file,append=T)
        }
    }
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

main(data_file,sep=sep,design_file=design_file, contrast_file=contrast_file, output_file=output_file)
