#!/usr/bin/env Rscript

# Usage: pca.R <input file> <output prefix> <columns>
# * input file: data table, one gene per row, one column per sample; first column contains gene ids.
# * output prefix: the result will be a pdf named "<prefix>.pdf"
# * columns:
#   - "all": selects all columns of the input table
#   - "rpkm": selects all columns which name contains "rpkm"
#   - "1,3,5,6...": comma-separated list of column indices (1-based, exclude the first column)
#
# The input table needs a header (column names).
# Column names are arbitrary, but are expected to be formatted as
# "[rpkm.]<group_name>.<#replicate>", e.g. "rpkm.Control.2" or "TRF2_KO.1".
# All samples with a different group name will get a different color.
#
# Ex: ./pca.R genes_expression.txt pca_biplot rpkm

filename = "genes_expression.txt"
outprefix = "pca_biplot"
columns = "all"

args = commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) { filename = args[1] }
if (length(args) >= 2) { outprefix = args[2] }
if (length(args) >= 3) { columns = args[3] }

# Columns selection
columns = tolower(columns)
if (columns == "all") {
    cols <- function(d) {1:ncol(d)}
} else if (columns == "rpkm") {
    cols <- function(d) {grep("rpkm.", colnames(d))}
} else {
    cols <- function(d) {unlist(lapply(strsplit(columns, split=','), FUN=function(x){as.numeric(x)}))}
}

d = read.table(filename, header=TRUE, row.names=1)
X = d[, cols(d)]
M = as.matrix(log(X+0.5))
#gene_names = d[,'GeneName']  # for arrow, TODO(?)
pca = prcomp(t(M), retx=TRUE, center=TRUE, scale.=FALSE)

# Limits
pcs = pca$x[,1:2]
xmin = min(pcs[,1])
xmax = max(pcs[,1])
ymin = min(pcs[,2])
ymax = max(pcs[,2])
margin = 0.4
xmin = xmin - abs(xmin*margin)
ymin = ymin - abs(ymin*margin)
xmax = xmax + abs(xmax*margin)
ymax = ymax + abs(ymax*margin)

# Group names
rnames = names(X)
rnames = sub("rpkm[.]", "", rnames)
gnames = sub("[.][0-9]$", "", rnames)

# Biplot and bar plot of variances
pdf(paste(outprefix,".pdf",sep=''), height=12, width=6)
par(mfrow=c(2,1), oma=c(0,0,0,0), pty="s", cex=0.75)
colors = as.numeric((as.factor(gnames)))
plot(pcs, pch=20, col=colors, xlim=c(xmin,xmax), ylim=c(ymin,ymax))
text(pcs, rnames, pos=3, col=colors)

vars = pca$sdev ^ 2
total_var = sum(vars)
barplot(100*(vars[1:10])/total_var, col=c(rep("grey",2),rep("white",8)),
    names.arg=sapply(1:10, FUN=function(x) paste("PC",x,sep='')),
    main=paste("PC1 + PC2: ", round(100*(sum(vars[1:2]))/total_var,2), "% of total variance")
)
dev.off()
