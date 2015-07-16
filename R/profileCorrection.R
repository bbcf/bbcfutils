
#profileCorrection.R --args infile baitCoord name outputTrack reportFile tableFile script_path
#inputFile, baitCoord, name, outputFile, reportFile, script_path=''

Args <- commandArgs(TRUE)
print(length(Args))
print("profileCorrection.R --args ")
print(Args)
infile <- Args[1]
baitcoord <- Args[2]
curName <- Args[3]
correctedFile <- Args[4]
report <- Args[5]
tableFile <- Args[6]

bin.log=function(x) {
    r1=min(x$start)
    r2=max(x$end)
    c((r1+r2)/2,sum(x$val*(x$end-x$start+1))/(1+r2-r1))
}
log.na=function(x,b=2){
    y=rep(NA,length(x))
    I=which(x>0)
    y[I]=log(x[I],base=b)
    y
}
unfold=function(x){
    cbind(unlist(apply(x,1,function(y)y[1]:y[2])),
          rep(x[,3],x[,2]-x[,1]+1))
}


profileCorrection <- function(fragsFile,baitCoord,plotFile=NA,all=FALSE)
{

        baitCoordSplit <- unlist(strsplit(gsub("-",":",baitCoord,perl=TRUE),":"))
        bait=as.numeric(baitCoordSplit[2:3])
        baitChr=baitCoordSplit[1]
        print(paste("Read file ",fragsFile,sep=""))
        bed=read.delim(fragsFile, skip=1,header=F,stringsAsFactors=F)
        if(ncol(bed)<4){bed=read.delim(fragsFile, skip=1,header=F,stringsAsFactors=F,sep=" ")}
        if(nrow(bed)==0){print("Warning!! File is empty!!")}

        I.down=which(bed[,1]==baitChr & bed[,2]>bait[2])
        data.down=data.frame(start=bed[I.down,2]-bait[2],end=bed[I.down,3]-bait[2],val=bed[I.down,4])
        I.up=which(bed[,1]==baitChr & bed[,3]<bait[1])
        data.up=data.frame(start=bait[1]-bed[I.up,3],end=bait[1]-bed[I.up,2],val=bed[I.up,4])
        cuts=round(log.na(data.down$start+data.down$end))
        binned.down=as.data.frame(t(sapply(split(data.down[which(data.down[,3]>0.0),],cuts[which(data.down[,3]>0.0)]),bin.log)))
        (lm.down=lm(y~x,data=data.frame(x=log(binned.down[,1]),y=log(binned.down[,2]))))
        cuts=round(log.na(data.up$start+data.up$end))
        binned.up=as.data.frame(t(sapply(split(data.up[which(data.up[,3]>0.0),],cuts[which(data.up[,3]>0.0)]),bin.log)))
        (lm.up=lm(y~x,data=data.frame(x=log(binned.up[,1]),y=log(binned.up[,2]))))
        a.down=mean(log(binned.down[,2]))+mean(log(binned.down[,1]))
        a.up=mean(log(binned.up[,2]))+mean(log(binned.up[,1]))
        cs.down<-1+exp(a.down)/ (0.5*(bed[I.down,2]+bed[I.down,3])-bait[2])
        d.down<-data.frame("start"=bed[I.down,2],"end"=bed[I.down,3],"score"=bed[I.down,4],"fit"=cs.down,"scoreFit"=bed[I.down,4]/cs.down)
        cs.up<-1+exp(a.up)/ (bait[1]-0.5*(bed[I.up,2]+bed[I.up,3]))
        d.up<-data.frame("start"=bed[I.up,2],"end"=bed[I.up,3],"score"=bed[I.up,4],"fit"=cs.up,"scoreFit"=bed[I.up,4]/cs.up)
        data2=rbind(d.up,d.down)
        data2=cbind(rep(baitChr,nrow(data2)),data2)
        colnames(data2)[1]="chromosome"

        if(!is.na(plotFile))
        {
                pdf(file=plotFile)
                plot(data2$start,data2$score,t='h',ylim=c(-max(data2$scoreFit),max(data2$score)),main=paste("corrected profiles\n",plotFile,sep=""),xlab="",ylab="frags score",col="dark grey")
                points(data2$start,-1*data2$scoreFit,t='h',col="orange")
                lines(0.5*(data2$start+data2$end),data2$fit,col="red",lty=2) #orangered or orangered1
                legend("topright",legend=c("before correction","profile corrected","fit"),col=c("darkgrey","orangered","red"),lwd=c(1,1,2))
                dev.off()


                #pdf(file=gsub(".pdf","_fit.pdf",plotFile))
                pdf(file=plotFile)
                plot(binned.down[,1],binned.down[,2],log='xy',pch=20,xlab='position',ylab='value',t='b',col='red')
                a.down=mean(log(binned.down[,2]))+mean(log(binned.down[,1]))
                lines(binned.down[,1],1+exp(a.down)/binned.down[,1],col='magenta')
                points(binned.up[,1],binned.up[,2],pch=20,col='blue',t='b')
                a.up=mean(log(binned.up[,2]))+mean(log(binned.up[,1]))
                lines(binned.up[,1],1+exp(a.up)/binned.up[,1],col='cyan')
                legend("topright",legend=c("Down (data)","Up (data)","Down (fit)","Up (fit)"),col=c("red","blue","magenta","cyan"),lty=c(2,2,1,1))
                title("data vs. fit")

                dev.off()
        }

    if(all)
    {   # will return data for all chromosomes
        I.otherChr <- which(bed[,1] != baitChr)
        res = cbind(bed[I.otherChr,],rep(NA,length(I.otherChr)),bed[I.otherChr,4])
        #colnames(res) = colnames(data2)[c(1:3,6)]
        colnames(res) = colnames(data2)
        #res <- rbind(res,data2[,c(1:3,6)])
        res <- rbind(res,data2)
        return(res)
}

    #return(data2[,c(1:3,6)]) #return only profile corrected data (=data on baitChr only)
    return(data2) #return only profile corrected data (=data on baitChr only)
}

correctedData=profileCorrection(infile,baitcoord,report,TRUE)
# was c(1:3,6)

header=paste("track type=bedGraph name='",curName," (profile corrected)' description='",curName," (profile corrected)' visibility=full windowingFunction=maximum",sep="")
write.table(header,correctedFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(cbind(correctedData[,1:3],round(correctedData[,6],2)),file=correctedFile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)

## export full table (score, fit, score/fit
write.table(cbind(correctedData[,1:3],round(correctedData[,4:6],2)),file=tableFile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

print("*****************")
print(paste("Profile Correction of ",infile," done!"))
print(paste("resfile=",correctedFile))
print(paste("tablefile=",tableFile))
print(paste("report file=",report))
print("*****************")


