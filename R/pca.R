#!/usr/bin/env Rscript

# Usage: pca.R <input file> <output prefix>


args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    filename = args[1]
    if (length(args) > 1) {
        outprefix = args[2]
    } else {
        outprefix = "pca_biplot"
    }
} else {
    filename = "genes_expression.txt"
}

d = read.table(filename, header=TRUE, row.names=1)
X = d[, grep("rpkm.", colnames(d))]
X = as.matrix(log(X+0.5))
gene_names = d[,'GeneName']
pca = prcomp(t(X), retx=TRUE, center=TRUE, scale.=FALSE)

# Limits
pcs = pca$x[,1:2]
xmin = min(pcs[,1])
xmax = max(pcs[,1])
ymin = min(pcs[,2])
ymax = max(pcs[,2])
margin = 0.2
xmin = xmin - abs(xmin*margin)
ymin = ymin - abs(ymin*margin)
xmax = ymax + abs(ymax*margin)
ymax = ymax + abs(ymax*margin)

# Group names
rnames = names(d)[grep("rpkm.",names(d))]
rnames = sub("rpkm.", "", rnames)
gnames = sub(".[0-9]$", "", rnames)

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
