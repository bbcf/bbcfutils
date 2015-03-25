#!/usr/bin/env Rscript

# Runs limma Voom to get DE transcripts
# Usage: limma.R <filename> -s <separator> -o <output_prefix>

args = commandArgs(trailingOnly = TRUE)

data_file = args[1]
sep = args[grep("-s",args)+1]
output_prefix = args[grep("-o",args)+1]
# As an arg in bash, a tab should be passed as $'\t'. '\t' is translated to "\\t".
if (sep=="\\t") {sep='\t'}

d = read.table(data_file, header=TRUE, row.names=1, sep=sep, quote="", check.names=F)
cc = grep("^[Cc]ount[s]*[.]", colnames(d))
snames = colnames(d)[cc]
conds = as.factor(sapply(strsplit(snames,'.',fixed=T),function(x){l=length(x);paste(x[2:(l-1)],collapse='.')}))
groups = as.vector(unique(conds))

# Statistics need replicates to work
couples = combn(groups,2)
ctable = table(conds)
exclude = c()
for (i in 1:dim(couples)[2]){
    if ((ctable[[couples[1,i]]]==1) && (ctable[[couples[2,i]]]==1)) {
        exclude = c(exclude, i)
    }
}
if (! is.null(exclude)) {
    couples = couples[,-exclude]
}

# Design
design = model.matrix(~ 0+conds)
colnames(design) = sub("conds", "", colnames(design))
rownames(design) = snames

library(limma)

# Contrasts
comps = c()
for (i in 1:dim(couples)[2]){
    comps = c(comps, paste(couples[1,i],"-",couples[2,i], sep=''))
}
args = as.list(comps)
args$levels = design
names(args) = c(comps,"levels")
contrasts = do.call(makeContrasts, args)

# Counts matrix
countdata = d[,cc]
rownames(countdata) = rownames(d)

# Remove lines with all counts < 50
# "The limma-voom method assumes that rows with zero or very low counts have been removed."
# (from their userguide.pdf, 15.3 Differential expression, p.69)
cd = countdata[apply(countdata, 1, function(row) any(row > 50)), ]

# Normalization
estimateSizeFactorsForMatrix <- function(counts)  {
   loggeomeans <- rowMeans( log(counts) )
   apply( counts, 2, function(cnts)
      exp( median( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}
norm.factors = estimateSizeFactorsForMatrix(cd)
norm.counts = t(t(cd)/norm.factors)

# DE analysis
v <- voom(norm.counts, design)
fit <- lmFit(v,design)
fit2 = contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
for (coef in colnames(contrasts)) {
    tt = topTable(fit2, coef=coef, adjust.method="BH", number=nrow(v), sort.by="P")
    # Reannotate
    res = cbind(tt, d[tt$ID, c("Chrom","Start","End","Strand","GeneName")])
    comp = gsub(" ","", coef , fixed=TRUE)
    outname = paste(output_prefix,"_",comp,".txt", sep='')
    write(coef, outname)
    write.table(res, outname, append=TRUE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
}

