.libPaths(new="/home/epfl/bbcf/lib/R/x86_64")
source("/sfs1/home/frt/jrougemo/Data/BBCFPeaks/Peaks/deconv_fcts.R")
filename="peaks_allZTvsInput_merged.bed"
regs=read.delim(filename,header=F)
regions=split(regs[,-1],regs[,1])
chr="_CHR_"
rm(regs)

wigf=paste("Pooled_",chr,"_deconv.wig",sep='')
bedf=paste("Pooled_",chr,"_deconv.bed",sep='')
pdff=paste("Pooled_",chr,"_deconv.pdf",sep='')

filename="Pooled_macs_ZTvsTI_regions.bed.gz"
data = read.data(filename,c(chr),regions,offset=TRUE,gff=FALSE)
d = data.to.strands(c(chr),data)
rm(data)
counts = regions2counts(regions,d,threshold=-1)
pdf(file=pdff,title=chr,paper="a4",width=8,height=11)
if (length(counts[[1]][[chr]]) > 0) {
    plot.counts(counts[[1]][[chr]],d$prod[[chr]])
    nfact = 2*sum(d$prod[[chr]][,3])/sum(c(d$plus[[chr]][,3],d$minus[[chr]][,3]))
    lines(d$minus[[chr]][,1]*1e-3,cumsum(d$minus[[chr]][,3])*nfact,col='red')
    lines(d$plus[[chr]][,1]*1e-3,cumsum(d$plus[[chr]][,3])*nfact,col='blue')
}
ccf = cross.correlate(counts,threshold=0)
plot(ccf$lag,ccf$acf,t='l',xlab='Lag',ylab='Cross-correlation',main='Strand cross-correlation',ylim=c(0,1))
sol = inverse.solve(counts,mu=80,lambda=150,len=38,regul=1e-5,optimize=FALSE)
print(sol$par)
header=c("track type=bedGraph visibility=2 name=pooled_deconv",chr)
append=FALSE
write.wig(cbind(counts[[1]][[chr]][[1]]$pos,counts[[1]][[chr]][[1]]$pos,sol$sol[[chr]][[1]]$prob),header,wigf,offset=TRUE,append=append,cutoff=1e-3)
if (length(counts[[1]][[chr]])>1) {
    for (n in 2:length(counts[[1]][[chr]])) {
        write.wig(cbind(counts[[1]][[chr]][[n]]$pos,counts[[1]][[chr]][[n]]$pos,sol$sol[[chr]][[n]]$prob),header,wigf,offset=TRUE,append=TRUE,cutoff=1e-3)
    }
}
write.bed.sol(sol$sol[[chr]],counts[[1]][[chr]],header,bedf,append=append,offset=TRUE,gff=FALSE,cutoff=1e-3)
par(mfrow=c(4,2))
for (n in 1:length(counts[[1]][[chr]])) {
    if (sol$sol[[chr]][[n]]$value>.65) next
    plot.sol(counts[[1]][[chr]][[n]],sol$sol[[chr]][[n]],sol$par)
    title(sub=chr)
}
dev.off()
