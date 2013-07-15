# This script uses DESeq to look for differential expression in RNA-seq data
# of genomic features in different conditions.

# --------
#  Usage:
# --------
# R --slave -f negbin.test.R --args 'data_file' -s 'tab' -o 'output_file'

# *** Personal example:
# (On Vital-IT, .bashrc contains: alias R="bsub -qbbcf -Ip -XF R --vanilla")

# R --slave -f negbin.test.R --args 'tests/testing_files/genes_expression.tab' -s 'tab' -o 'tests/testing_files/output_deseq.txt'

# ------------
#  Arguments:
# ------------
# -s sep: character separating fields in the input files. Use 'tab' for tab delimiter (not '\t').
# -o output_file: name of the file(s) containing the results. As many files as there are pairs of groups are
# created - the corresponding pair is added as a suffix to the file name. If not specified, the result is
# printed to stdout.

# --------
#  Input:
# --------
# The script is made to take as input the data files returned by rnaseq.py, which format is the following:
# - tab-delimited or CSV file containing raw read counts
# - lines representing genomic features
# - the first column contains feature names
# - other columns contain counts, one column per sample
# - columns containing counts are labeled "counts.<groupName>.<sampleIndex>" (e.g. counts.G1.3).

library(MASS)

args = commandArgs(trailingOnly = TRUE)

data_file = args[1]
sep = args[grep("-s",args)+1]
output_file = args[grep("-o",args)+1]
if (sep=='tab') sep='\t'

#options(error = quote({dump.frames(to.file=TRUE); q()})) # creates an error log file `last.dump.rda`

main <- function(data_file, sep="\t", output_file=''){
    data = read.table(data_file, header=T, row.names=1, sep=sep, quote="", check.names=F)
    header = colnames(data)
    counts = grep("^counts[.]",header)
    data = round(data[,counts])
    samples = header[counts]
    conds = sapply(strsplit(samples,'.',fixed=T),function(x){l=length(x);paste(x[2:(l-1)],collapse='.')})

    # Still need to check that replicates are not identical - lfproc would fail
    if (all(table(conds)>3)){        # if >3 replicates in all conditions
        method = 'per-condition'        # for each group estimate the variance from its replicates
        sharingMode = 'gene-est-only'   # use the per-gene variance estimates only
    } else if (any(table(conds)>1)){ # if few replicates
        method = 'pooled'               # use all groups with replicates to estimate the variance
        sharingMode = 'maximum'         # use the max of the GLM fit and the estimated variance
    } else {                         # if no replicates
        method = 'blind'                # pools all groups together to estimate the variance
        sharingMode='fit-only'          # use only the GLM fit across the pooled variance
    }

    if (nrow(data)>3){
    DES(data, conds, method, sharingMode, output_file) }
}


# Run DESeq on *data*, the counts matrix
# Tests all possible 2-by-2 combinations of groups
DES <- function(data, conds, method, sharingMode, output_file=FALSE){
    library(DESeq)
    groups = unique(conds)

    result = list()
    cds <- newCountDataSet(data, conds)
    cds <- estimateSizeFactors(cds)
    test = try({
        cds <- estimateVarianceFunctions(cds, method=method)
    }, silent=TRUE)
    if(class(test) == "try-error") {
        test2 = try({
            cds <- estimateDispersions(cds, method=method, fitType='parametric', sharingMode=sharingMode)
        })
        if(class(test2) == "try-error") {
            cds <- estimateDispersions(cds, method=method, fitType='local', sharingMode=sharingMode)
        }
    }
    couples = combn(groups,2)
    for (i in 1:dim(couples)[2]){
        res <- nbinomTest(cds, couples[1,i], couples[2,i])
        res = res[order(res[,'padj']),] # sort w.r.t. adjusted p-value
        comp = paste(couples[1,i],"-",couples[2,i])
        result[[comp]] = res
    }
    write_result(output_file, result)
}


# Write an output table for each of the comparisons made
write_result <- function(output_file, res_list, sep='\t'){
    if (length(output_file) == 0){
        print(res_list)
    }else{
        for (x in names(res_list)){
            xs = unlist(strsplit(x,' '))
            xs = paste(xs[1],xs[2],xs[3],sep='')
            out = paste(output_file,xs,sep='_')
            write(x,out)
            write.table(res_list[[x]],out,quote=F,row.names=T,col.names=T,append=T,sep="\t")
        }
    }
}


main(data_file,sep=sep,output_file)


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#


