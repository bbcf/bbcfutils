options(stringsAsFactors=F)
library(rjson)
args=commandArgs(trailingOnly = TRUE)
stats.by.sample=fromJSON(file=args[1])
pdf(file=args[2],paper="a4",height=11,width=8)
par(cex=1.5,lwd=1.5,mfrow=c(4,1),oma=c(4,0,4,0))
for (sample in sort(names(stats.by.sample))) {
    stats=stats.by.sample[[sample]]
    df=data.frame(hits=names(stats$multi_hits),reads=as.numeric(stats$multi_hits))
    df=df[order(as.numeric(df$hits)),]
    if (stats$unmapped>0) df=rbind(df,c(0,as.numeric(stats$unmapped)))
#    col='darkorange'
    col=heat.colors(6)[c(2,rep(4,length(stats$multi_hits)-1),6)]
    read.total = stats$total+stats$unmapped
    p=barplot(df$reads+.1,names.arg=df$hits,border=0,
      log='y',col=col,xlab='# hits',ylab='# reads',ylim=c(1,read.total*10),
      main='Reads with multiple hits')
    text(x=p,y=5,lab=df$reads,srt=90,adj=c(0,0.5),cex=1.1)
    abline(h=stats$total)
    text(nrow(df)/3,stats$total,paste("reads mapped:",stats$total),pos=3)
    if (read.total>stats$total) {
        abline(h=read.total,lty=2)
        text(2*nrow(df)/3,read.total,paste("reads total:",read.total),pos=3)
    }
    df=data.frame(mismatches=names(stats$mismatches),
      reads=as.numeric(stats$mismatches))
    df=df[order(as.numeric(df$mismatches),decreasing=T),]
    if (nrow(df) > 9) df=rbind(data.frame(mismatches=">8",reads=sum(df[1:(nrow(df)-9),"reads"])),
                        df[nrow(df)-8:0,]) 
    col=heat.colors(6)[c(rep(4,length(df$reads)-1),2)]
#    col='darkorange'
    p=barplot(df$reads,names.arg=df$mismatches,horiz=T,border=0,
      col=col,ylab='# mismatches',xlab='# reads',las=1,
      main='Mismatches between reads and reference')
    text(x=mean(df$reads),y=p,lab=df$reads,cex=1.1)
    ps=as.numeric(c(stats$alignments$fwd,stats$alignments$rev))
    if (sum(ps)>0) {
      col=c('blue','red')
      pie(ps,labels=paste(c("Forward","Reverse"),ps),col=col,
          main='Alignements per strand of reference',cex=1.1)
    }
    df=data.frame(genome=as.numeric(stats$genome_size),
      expected=as.numeric(stats$expected_coverage),
      actual=as.numeric(stats$actual_coverage))
    col=c('darkorange','forestgreen','forestgreen')
    p=barplot(c(df$genome,df$genome*df$expected,df$genome*df$actual)*1e-6,
      horiz=T,border=0,col=col,ylab='',xlab='Megabase',main='Genome coverage')
    text(x=df$genome*5e-7,y=p,lab=c(
                   sprintf("genome size: %.0f",df$genome),
                   sprintf("total read length: %.1f%%",df$expected*100),
                   sprintf("actual portion covered: %.1f%%",df$actual*100)),cex=1.1)
    title(main=sample,outer=T,cex.main=4)
    if (length(stats$cmd_line)) mtext(stats$cmd_line,1,outer=T,cex=.8,adj=0)
}
dev.off()
