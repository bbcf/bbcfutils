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

args=commandArgs(trailingOnly = TRUE)

data_file = args[1]
sep = args[grep("-s",args)+1]
output_file = args[grep("-o",args)+1]
if (sep=='tab') sep='\t'

#options(error = quote({dump.frames(to.file=TRUE); q()})) # creates an error log file `last.dump.rda`

main <- function(data_file, sep="\t", output_file=''){
    data = read.table(data_file, header=T, row.names=1, sep=sep)
    header = colnames(data)
    counts = grep("counts",header,fixed=T)
    data = round(data[,counts])
    samples = header[counts]
    conds = unlist(lapply(strsplit(samples,".",fixed=T), "[[", 2))

    # Still need to check that replicates are not identical - lfproc would fail
    if (any(table(conds)>1)){ method = 'normal' # if replicates
    } else { method = 'blind' }

    if (nrow(data)>3){
    DES(data, conds, method, output_file) }
}

write_result <- function(output_file, res_list, sep='\t'){
    if (length(output_file) == 0){
        print(res_list)
    }else{
        for (x in names(res_list)){
            xs = unlist(strsplit(x,' '))
            xs = paste(xs[1],xs[2],xs[3],sep='')
            output_file = paste(output_file,xs,sep='_')
            write(x,output_file)
            write.table(res_list[[x]],output_file,quote=F,row.names=T,col.names=T,append=T,sep="\t") # write.csv?
        }
    }
}

DES <- function(data, conds, method='normal', output_file=FALSE){  ## DESeq ##
    library(DESeq)
    groups = unique(conds)

    result = list()
    cds <- newCountDataSet(data, conds)
    cds <- estimateSizeFactors(cds)
    cds <- estimateVarianceFunctions(cds, method=method)
    couples = combn(groups,2)
    contrast.names = c()
    for (i in 1:dim(couples)[2]){
        res <- nbinomTest(cds, couples[1,i], couples[2,i])
        res = res[order(res[,8]),] # sort w.r.t. adjusted p-value
        comp = paste(couples[1,i],"-",couples[2,i])
        contrast.names = c(contrast.names,comp)
        result[[comp]] = res
    }
    write_result(output_file, result)
}


main(data_file,sep=sep,output_file)


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

