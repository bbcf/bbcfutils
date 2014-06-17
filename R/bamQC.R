#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
pdfname = args[1]
gsize = args[2]
title = args[3]

data = read.delim(file('stdin'), header=FALSE)
freqs = sapply(split(data[,3]-data[,2],data[,4]),sum)
freqs0 = gsize-sum(freqs)
qtiles = freqs0+cumsum(freqs)
med = as.integer(names(freqs)[max(which(qtiles < gsize/2 ))])
avg = round(sum(as.numeric(names(freqs))*freqs)/gsize,2)

pdf(pdfname)
par(cex=1.3,lwd=1.5,las=1,mar=c(5,5,5,1)+.1,mgp=c(4,1,0))
plot(as.integer(names(freqs)),freqs,log='xy',pch=20,t='b',ylim=c(1,gsize*1.5),
     xlab='Nb reads',ylab='Nb genomic positions',main=title)
abline(h=c(freqs0,gsize),col=c('blue','red'))
text(2,freqs0,paste("Uncovered positions:",freqs0),pos=1,col='blue',cex=.8)
text(2,gsize,paste("Genome size:",gsize),pos=3,col='red',cex=.8)
abline(v=c(avg,med),col='cyan',lty=2)
text(med,1,paste("Median:",med),pos=2,col='cyan',cex=.8)
text(avg,5,paste("Average:",avg),pos=2,col='cyan',cex=.8)
dev.off()

