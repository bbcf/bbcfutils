args = commandArgs(trailingOnly = TRUE)
counts.file = args[1]
pdf.file = args[2]
read.length = as.integer(args[3])
chr.name = args[4]
output.file = args[5]
script.path = args[6]
mu = 80
options(stringsAsFactors=F)
source(paste(script.path,"/deconv_fcts.R",sep=''))

counts = read.data(counts.file)
pdf(file=pdf.file,title='chip-seq',paper="a4",width=8,height=11)
par(cex=1.5,lwd=1.5)
ccf = cross.correlate(counts,threshold=.5)
plot(ccf$lag,ccf$acf,t='l',ylim=c(0,1),
     xlab='Lag',ylab='Cross-correlation',
     main=paste('Strand cross-correlation',chr.name))
cut.ccf=ccf$acf
cut.ccf[which(ccf$lag<mu)]=0
lambda=ccf$lag[which.max(cut.ccf)]
sol = inverse.solve(counts,mu=mu,lambda=lambda,len=read.length,regul=1e-3,optimize=TRUE)
col = 'red'
#lab=substitute(expression(lambda=x),list(x=sol$par$lambda))
lab = paste('lambda=',sol$par$lambda,sep='')
abline(v=sol$par$lambda,col=col)
text(sol$par$lambda,0,lab,col=col,pos=4)
col = 'blue'
#lab=substitute(expression(mu=x),list(x=sol$par$mu))
lab = paste('mu=',sol$par$mu,sep='')
abline(v=sol$par$mu,col=col)
text(sol$par$mu,0.3,lab,col=col,pos=4)
col = 'darkgreen'
#lab=substitute(expression(ell=x),list(x=read.length))
lab = paste('l=',read.length,sep='')
abline(v=read.length,col=col)
text(read.length,0.6,lab,col=col,pos=4)
par(mfrow=c(4,2))
for (n in names(counts)) {
    if (sol$sol[[n]]$value>.65) next
    plot.sol(counts[[n]],sol$sol[[n]],sol$par)
    title(sub=chr.name)
}
dev.off()
bed = data.frame()
cutoff = 1e-3
for (n in names(counts)) {
    I = which(sol$sol[[n]]$prob>cutoff*sum(sol$sol[[n]]$prob))
    if (length(I)<2) next
    interval = range(counts[[n]]$pos[I])
    score = sum(sol$sol[[n]]$prob[I])
    name = paste("ID=",n,";FERR=",round(sol$sol[[n]]$val,digits=4),sep='')
    bed = rbind(bed,
      data.frame(chr=chr.name,
                 start=interval[1],
                 end=interval[2],
                 name=name,
                 score=score))
}
bed[,"start"] = bed[,"start"]-1
wig = data.frame()
for (n in names(counts)) {
    I = which(sol$sol[[n]]$prob>cutoff*sum(sol$sol[[n]]$prob))
    wig = rbind(wig,data.frame(
      chr=rep(chr.name,length(I)),
      pos=as.integer(counts[[n]]$pos[I]),
      score=as.numeric(sol$sol[[n]]$prob[I])))
}
par = sol$par
save(bed,wig,par,file=output.file)
