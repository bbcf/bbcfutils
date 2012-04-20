
# identify interactive region: approach 2
predictInteractiveRegion <- function(myData, baitCoord, scale=2000000,myPercentile=0.95,myTitle="",plotFile)
{
	i_bait_upstream <- which(myData[,1]==baitCoord[,1] & myData[,2] > as.numeric(baitCoord[,3]))
	i_bait_downstream <- which(myData[,1]==baitCoord[,1] & myData[,3] < as.numeric(baitCoord[,2]))
	bait_center <- as.numeric(baitCoord[,2])+((as.numeric(baitCoord[,3])-as.numeric(baitCoord[,2]))/2)

	############################		
	### Do it for Downstream region
	############################		
	x <- abs(myData[i_bait_downstream,2]-bait_center)
	y <- myData[i_bait_downstream,4]	
	#i=1:length(x) 
	i <- which(x < scale)

	pdf(file=plotFile)	
	par(mfrow=c(2,1))
	
	#partI: plot original fit of log (verify that eq=a+b*d^(-2/3) can be applied.)
	#fitLog$coefficients[2] should be close to -2/3
	plot(log(x[i]),log(y[i]),pch=20,cex=0.8,xlab="distance to the center of the bait (log)",ylab="score (log)",main=paste(myTitle,"\nSelection of Interacting Region (downstream)"))
	fitLog<- lm(a~b,data=list(a=log(0.1+y[i]),b=log(0.1+x[i])))
	lines(log(0.1+x[i]),predict(fitLog),col="blue")
	text(x=min(log(0.1+x[i]))+1,y=min(predict(fitLog)),paste("coeff=",round(fitLog$coefficients[2],2)),col="blue")
	print(fitLog$coefficients[2])

	#partII: deduce dmax such that eq=a+b*d^(-2/3) < th
	#th set as 90% percentile of pred(fit), "fit" being the linear model y=a+b*d^(-2/3)
	coeff=-1 #-2/3 #-1
	fit <- lm(a~b+0,data=list(a=y[i],b=(1+x[i])^coeff)) # a~b+0 => no intercept
	z <- predict(fit)
	
	boxplot(z,ylim=c(0,2000))
	plot(density(z))
	##th<-0.1*max(z)
#	th <- quantile(z,myPercentile)
	bckgReg=scale/2
	j <- union(which(myData[,3] < (bait_center-(bckgReg/2))),which(myData[,2] > (bait_center+(bckgReg/2))) )
	j2 <- which(myData[j,2]> (bait_center-(3*scale)) & myData[j,3]< (bait_center+(3*scale)))
	th <- quantile(myData[j[j2],4],myPercentile)
	
	#th <- 2*IQR(z)
	invCoeff=-1 #-3/2 #-1
	#dmax <- ((th-fit$coefficients[1])/fit$coefficients[2])^invCoeff
	dmax <- (th/fit$coefficients)^invCoeff
	print(baitCoord)
	interactiveReg <- c(as.character(baitCoord[,1]),round(as.numeric(baitCoord[,3]),0),round(as.numeric(baitCoord[,3])+dmax,0))
	print(paste("th=",th,"=> dmax=",dmax))
	print(paste("interactive region downstream=",interactiveReg))
	
	plot(log(0.1+x[i]),z,type="b",pch=1,xlab="log(distance to bait)",ylab="predicted scores")
	abline(h=th,col="red",lty=2)
	text(min(log(0.+x[i]))+0.5,y=th+600,paste("th=",round(th,0)),col="red")
	abline(v=log(dmax),col="red",lty=2)
	text(x=log(dmax)+1,y=(max(z)-min(z))/2,paste("dmax=",round(dmax,0),"bps (log=",round(log(dmax),1),")"),col="red")
title("predicted scores vs. distance to bait")

	res <- interactiveReg

	############################		
       ### Do it for Upstream region
	###########################
        x <- abs(myData[i_bait_upstream,2]-bait_center)
        y <- myData[i_bait_upstream,4]
        i=1:length(x) #i <- which(x < scale)

        par(mfrow=c(2,1))

        #partI: plot original fit of log (verify that eq=a+b*d^(-2/3) can beapplied.)
        #fitLog$coefficients[2] should be close to -2/3
        plot(log(x[i]),log(y[i]),pch=20,cex=0.8,xlab="distance to the center of the bait (log)",ylab="score (log)",main=paste(myTitle,"\nSelection of Interacting Region (Upstream)"))
        fitLog<- lm(a~b,data=list(a=log(0.1+y[i]),b=log(0.1+x[i])))
        lines(log(0.1+x[i]),predict(fitLog),col="blue")
        text(x=min(log(0.1+x[i]))+1,y=min(predict(fitLog)),paste("coeff=",round(fitLog$coefficients[2],2)),col="blue")
        print(fitLog$coefficients[2])

        #partII: deduce dmax such that eq=a+b*d^(-2/3) < th
        #th set as 90% percentile of pred(fit), "fit" being the linear model y=a+b*d^(-2/3)
	coeff=-1 #-2/3 #-1	
        fit <- lm(a~b+0,data=list(a=y[i],b=(1+x[i])^coeff)) # a~b+0 => no intercept
        z <- predict(fit)
	boxplot(z,ylim=c(0,2000))
        plot(density(z))
	
	##th<-0.1*max(z)
#	th <- quantile(z,myPercentile)
	#th <- 2*IQR(z)
	invCoeff=-1 # -3/2 #-1
        #dmax <- ((th-fit$coefficients[1])/fit$coefficients[2])^invCoeff
        dmax <- (th/fit$coefficients)^invCoeff
        print(baitCoord)
        interactiveReg <- c(as.character(baitCoord[,1]),round(as.numeric(baitCoord[,2])-dmax,0),baitCoord[,2])
        print(paste("th=",th,"=> dmax=",dmax))
	print(paste("interactive region upstream=",interactiveReg))

        plot(log(0.1+x[i]),z,type="b",pch=1,xlab="log(distance to bait)",ylab="predicted scores")
        abline(h=th,col="red",lty=2)
        text(min(log(0.+x[i]))+0.5,y=th+600,paste("th=",round(th,0)),col="red")
        abline(v=log(dmax),col="red",lty=2)
	text(x=log(dmax)+1,y=(max(z)-min(z))/2,paste("dmax=",round(dmax,0),"bps (log=",round(log(dmax),1),")"),col="red")
	title("predicted scores vs. distance to bait")

        res[2] <- interactiveReg[2]
	interactiveReg <- res
	print(paste(" combined interactive region=",interactiveReg))


	########################
	# Plot Final Interactive Region
	#######################
	par(mfrow=c(2,1))
	
	my_ylim=5000
        my_xlim=c(as.numeric(baitCoord[,2])-2*scale,as.numeric(baitCoord[,3])+2*scale)
        plot(myData[,2],myData[,4],col="grey",type="h",xlim=my_xlim,ylim=c(0,my_ylim),ylab="fragment score",xlab="")
        print(interactiveReg)
        i_reg <- which(as.character(myData[,1])==as.character(interactiveReg[1]) & as.numeric(myData[,2])>=as.numeric(interactiveReg[2]) & as.numeric(myData[,3])<=as.numeric(interactiveReg[3]))
        i_bait <- which(myData[,1]==baitCoord[,1] & myData[,2]>=baitCoord[,2] & myData[,3]<=baitCoord[,3])
        points(myData[i_reg,2],myData[i_reg,4],col="orange",type="h")
        points(myData[i_bait,2],myData[i_bait,4],col="red",type="h")
        text(x=as.numeric(baitCoord[,3])+0.25*scale,y=2*my_ylim/3,paste("interacting region:\n",interactiveReg[1],":",interactiveReg[2],"-",interactiveReg[3],"\n(length=",(as.numeric(interactiveReg[3])-as.numeric(interactiveReg[2]))/1000,"kbps)"))
	plot.new()


        print("plot Cumulative scores (excluding self-ligated fragment)")
        plotCumSum_RegionFromCenter_exclSelfLigated(myData,baitCoord[,1],myCenter=round(as.numeric(bait_center),0),myExtension=0.5*scale,reg=interactiveReg,showReg=TRUE,myCol="blue",toExclude=baitCoord,add=FALSE,plotFile=plotFile,add2pdf=TRUE,titlePlot=paste("Cumulative distribution of scores\n +/- ",round(0.5*scale,0)/1000000,"Mbps around the center"))
        print("plot Cumulative scores for interacting region (excluding self-ligated fragment)")
        data_selfLigated <- plotCumSum_RegionFromCenter_exclSelfLigated(myData,baitCoord[,1],myCenter=round(as.numeric(bait_center),0),myExtension=dmax,myCol="blue",toExclude=baitCoord,add=FALSE,plotFile=plotFile,add2pdf=TRUE,titlePlot=paste("Cumulative distribution of scores in interacting region:\n",interactiveReg[1],":",interactiveReg[2],"-",interactiveReg[3]))

        dev.off()

	return(interactiveReg)
}

# identify interactive region: approach 1
plotMeanvsWinSize <- function(myData,baitCoord,step=5000,maxSize=500000,scale=1000000,myThType="IQR",myTitle="",plotFile)
{
	i_baitChr <- which(myData[,1] == baitCoord[,1])
	i_bait <- which(myData[,1]==baitCoord[,1] & myData[,2] >= as.numeric(baitCoord[,2]) & myData[,3] <= as.numeric(baitCoord[,3]))
	s=as.numeric(baitCoord[,2])
	e=as.numeric(baitCoord[,3])
	m=maxSize/(2*step)
	avgScore <- rep(0,m)
	winSize <- rep(0,m)
	k=1
	while(e-s < maxSize)
	{
		i <- which(myData[i_baitChr,2] >= s & myData[i_baitChr,2] <= e)
		avgScore[k] <- mean(myData[i_baitChr[i],4]) 
		winSize[k] <- e-s
		s <- s-step
		e <- e+step
		k=k+1
	}
	if(myThType=="IQR")
	{
	th=1.5*IQR(avgScore[-1],na.rm=TRUE)
	}
	else
	{	
	#th=2*IQR(avgScore[-1],na.rm=TRUE)
	th=quantile(avgScore[-1],as.numeric(myThType),na.rm=TRUE)
	}
	pdf(file=plotFile)	
	par(mfrow=c(2,1))
	plot(winSize[-1],avgScore[-1],cex=0.8,pch=16,xlab="region size (centered on bait)",ylab="mean score",main=paste(myTitle,"\nSelection of Interacting Region"))
	x <- which(avgScore[-1]>th)
	points(winSize[x],avgScore[x],col="red")
	abline(h=th,col="blue",lty=3)
	text(x=min(winSize[-1])+(2*(max(winSize[-1])-min(winSize[-1]))/3),y=th+1000,col="blue",sprintf("threshold=%.1f",th))
	
	s=as.numeric(baitCoord[,2])-(winSize[max(x)]/2)
	e=as.numeric(baitCoord[,3])+(winSize[max(x)]/2)
	
	my_ylim=5000
	plot(myData[i_baitChr,2],myData[i_baitChr,4],type="h",col="grey",ylim=c(0,my_ylim),xlim=c(as.numeric(baitCoord[,2])-scale,as.numeric(baitCoord[,3])+scale),ylab="fragment score",xlab="")
	x <- which(myData[i_baitChr,2]>s & myData[i_baitChr,3]< e)
	points(myData[i_baitChr[x],2],myData[i_baitChr[x],4],col="red",type="h")
	text(x=e+(scale/2),y=3*my_ylim/4,sprintf("interacting region\n%s:%s-%s\n (length=%.2f kbps)",baitCoord[,1],s,e,(e-s)/1000))
	dev.off()	

	return(c(as.character(baitCoord[,1]),as.numeric(s),as.numeric(e)))

}

#################################################
groupFrag <- function(myVector,step)
{
	s=seq(1,length(myVector)-step+1,by=step); 
	e=seq(step,length(myVector),by=step);
	return(lapply(1:length(s),function(j){mean(myVector[s[j]:e[j]]) }))
}

groupFrag_betweenIndices <- function(myVector,s,e)
{
        return(lapply(1:length(s),function(j){mean(myVector[s[j]:e[j]]) }))
}


generateWindows <- function(myData,n)
{
	data_byChr <- splitData(myData)
	#data_byChr <- split(myData,myData[,1])
	newWindows <- lapply(1:length(data_byChr), 
		function(i){
			s=c();e=c();
			if(nrow(data_byChr[[i]])>0){
			if(nrow(data_byChr[[i]]) < n){s=c(1);e=c(nrow(data_byChr[[i]]))}
			else{
			s=seq(1,nrow(data_byChr[[i]]),by=n); e=seq(n,nrow(data_byChr[[i]]),by=n);
#			s=seq(1,nrow(data_byChr[[i]]),by=1); e=seq(n,nrow(data_byChr[[i]]),by=1);
			if(length(e)<length(s)){e <- c(e,nrow(data_byChr[[i]]))}
			#if(length(e)<length(s)){e <- c(e,rep(nrow(data_byChr[[i]]),length(s)-length(e)))}

			}
			cbind(data_byChr[[i]][s,2],data_byChr[[i]][e,3],groupFrag_betweenIndices(data_byChr[[i]][,4],s,e))
#			cbind(data_byChr[[i]][s,2],data_byChr[[i]][e,3],groupFrag(data_byChr[[i]][,4],1))
			}
		}
	)
	names(newWindows)=names(data_byChr) 
	return(newWindows)
}

generateSmoothedWindows <- function(myData,n)
{
        data_byChr <- splitData(myData)
	if(length(data_byChr)>0){
        newWindows <- lapply(1:length(data_byChr),
                function(i){if(nrow(data_byChr[[i]])>0){
			s=c();e=c();
			mid=round(n/2);if((n/2)==round(n/2)){nprev=mid-1;nnext=mid}else{nprev=mid;nnext=mid}
			mid=seq(1,nrow(data_byChr[[i]]),by=1)
			istart=mid-rep(nprev,length(mid))
			iend=mid+rep(nnext,length(mid))
			ipairs <- cbind(istart,iend,mid)
			ipairs[which(ipairs[,2]>length(mid)),2]=length(mid)
			ipairs[which(ipairs[,1]<0),1]=1
                        
                        cbind(data_byChr[[i]][mid,2],data_byChr[[i]][mid,3],groupFrag_betweenIndices(data_byChr[[i]][,4],ipairs[,1],ipairs[,2]))
			} #end if nrow>0
                        } #end function(i) 
                
        )
        names(newWindows)=names(data_byChr)}
	else{newWindows=c()}
        return(newWindows)
}

# mid=seq(1,nrow(data_byChr[[i]]),by=1)
# n=18;mid=round(n/2);if((n/2)==round(n/2)){nprev=mid-1;nnext=mid}else{nprev=mid;nnext=mid}
#istart=mid-rep(nprev,length(mid))
#iend=mid+rep(nnext,length(mid))
#ipairs <- cbind(istart,iend,mid)
#> ipairs[which(ipairs[,2]>length(mid)),2]=length(mid)
#> ipairs[which(ipairs[,1]<0),1]=1

updateFactors <- function(myData)
{
	myData[] <- lapply(myData, function(x) x[,drop=TRUE])
	return(myData)
}

#getRegions <- function(myData,baitCoord,interRegCoord)
getRegions <- function(myData,interRegCoord)
{
	i_baitChr <- which(myData[,1] == interRegCoord[1])
	notBaitChr <- updateFactors(myData[which(myData[,1] != interRegCoord[1]),])
	i <- which(myData[i_baitChr,2] < as.numeric(interRegCoord[2]))
	regUp <- updateFactors(myData[i_baitChr[i],])
	i <- which(myData[i_baitChr,3] > as.numeric(interRegCoord[3]))
	regDown <- updateFactors(myData[i_baitChr[i],])
	return(list(othersChrs=notBaitChr,upstream_of_interactiveReg=regUp,down_of_interactiveReg=regDown))
}

mergeLists <- function(x,y)
{
	if (length(x) == 0){return(y)}
    if (length(y) == 0){return(x)}
	i = intersect(names(y), names(x))
	j <- setdiff(names(x),names(y))
	k <- setdiff(names(y),names(x))
	z <- list()
	if(length(i)>0){for(n in i){z[[n]] <- rbind(x[[n]],y[[n]])}}
	if(length(j)>0){for(n in j){z[[n]]<-x[[n]]}}
	if(length(k)>0){for(n in k){z[[n]]<-y[[n]]}}
	return(z)
}


splitData <- function(myData)
{
	z <- list();
	for(i in names(table(myData[,1]))){j <- which(myData[,1]==i);z[[i]]<-cbind(myData[j,],rownames(myData)[j])}
	return(z)
}

#----------------------------
plotCumSum_RegionFromCenter_exclSelfLigated <- function(myData,myChr,myCenter,myExtension,myCol,toExclude,reg="",showReg="FALSE",add="FALSE",plotFile="",plot="TRUE",add2pdf="FALSE",titlePlot="Cumulative distribution of scores")
{
myStart <- myCenter-myExtension
myEnd <- myCenter+myExtension
myIndex <- which(as.character(myData[,1]) == myChr & as.numeric(myData[,2])>=as.numeric(myStart) & as.numeric(myData[,3])<=as.numeric(myEnd))
myIndexSelfLigated <- intersect(which(as.numeric(myData[,1]) == as.character(toExclude[1]) & as.numeric(myData[,3])>=as.numeric(toExclude[2])),which(as.character(myData[,1])==as.character(myChr) & as.numeric(myData[,2])<=as.numeric(toExclude[3])))
myIndex <- setdiff(myIndex,myIndexSelfLigated)
print(paste("plotted:",length(myIndex)," (",length(myIndexSelfLigated)," frags excluded (self-ligated)",sep=""))
print(paste("Have excluded:",myIndexSelfLigated,sep=""))
if(plot==TRUE)
{
if(add2pdf==FALSE){pdf(file=plotFile)}
par(mfrow=c(1,1))
if(add==TRUE)
{
	points(myData[myIndex,2]-myCenter,cumsum(myData[myIndex,4]/sum(myData[myIndex,4])),cex=0.3,col=myCol,type="l")
}
else
{
	plot(myData[myIndex,2]-myCenter,cumsum(myData[myIndex,4]/sum(myData[myIndex,4])),cex=0.3,col=myCol,type="l",ylab="cumulative fragments scores",xlab="distance to center",main=titlePlot)
	abline(v=0,col="red",lty=3,lwd=1.2)
	if(showReg==TRUE)
	{
		abline(v=as.numeric(reg[2])-myCenter,col="orange",lty=2,lwd=1.2)
		abline(v=as.numeric(reg[3])-myCenter,col="orange",lty=2,lwd=1.2)
		legend("topleft",c(paste("center=",myChr,":",myCenter,sep=""),paste("interacting region\n",myChr,":",reg[2],"-",reg[3],sep="")),box.col="white",lty=c(3,2),lwd=c(1.2,1.2),col=c("red","orange"))
	}
	else{legend("topleft",paste("center=",myChr,":",myCenter,sep=""),box.col="white",lty=3,lwd=1.2,col="red")}
}
if(add2pdf==FALSE){dev.off()}
}
return(myIndexSelfLigated)
}


#--------------------------
findMinScore_FRD_perChr <- function(myData,myData_random,myFDR)
{
        print(paste("in findMinScore_FRD_perChr with FDR=",myFDR))
        #allScores <- c(myData,myData_random)
        allScores <- myData
        allScores_sorted <- sort(allScores,decreasing = TRUE)
        prevTh= -1
        for(i in 1:length(allScores_sorted))
        {
                curTh <- allScores_sorted[i]
                nOk <- length(which(myData>=curTh))
                nOkInRand <- length(which(myData_random>=curTh))
                a=0; #0.001
                curFDR <- ((nOkInRand+a)/length(myData_random))/((nOk+a)/length(myData))
                print(paste("curTh=",curTh,"nOkInRand=",(nOkInRand/length(myData_random)),"Ok=",(nOk/length(myData)),"=> curFDR=",curFDR))
                if(curFDR>myFDR){print(paste("th found=",prevTh,"=>",nOk,"satisfy this threshold")); return(c(prevTh,nOk,sprintf("%.2f",nOk/length(myData))))}
                prevTh=curTh
        }
        print("No th found")
        return(c("NA","0","0"))

}

findMinScore_FDR <- function(myData,nFragsPerWin,myNsamples,myFDR)
{
        allFoundTh <- lapply(1:length(myData),
        function(i)
        {
                print(paste("will treat chr",names(myData)[i]));
                print(class(unlist(myData[[i]][,3])))
                nWindows = length(unlist(myData[[i]][,3]))
                myNsamples=10*nWindows
                print(paste("will generate background composed by ",myNsamples," random windows (for",nWindows," windows on chr",names(myData)[i],"dim(myData)=",dim(myData[[i]])))
                if(length(unlist(myData[[i]][,3]))>=nFragsPerWin)
		{
		background <- unlist(lapply(1:myNsamples,function(j){mean(sample(unlist(myData[[i]][,3]),nFragsPerWin))}))
                findMinScore_FRD_perChr(unlist(myData[[i]][,3]),background,myFDR)
		}
		else
		{	
			print("#windows < nFragsPerWindows => no Th found")
			c(max(unlist(myData[[i]][,3])),0,0)	
		}
        })
        names(allFoundTh) = names(myData)
        return(allFoundTh)
}


#---------------------------
writeWindows <- function(myData,mySelectedTh,resFile,selectionFile)
{
        x <- splitData(mySelectedTh)
        curStatus <- lapply(1:length(x),
        function(i)
        {
		curChr=names(x)[i];	
                curTh=as.numeric(unlist(x[[i]][2])); print(paste(curChr,"=>",curTh)); 
		y=unlist(myData[[curChr]][,3]);
                i <- which( (y+0.001) >=curTh)
		curRes <- cbind(curChr,myData[[curChr]],(y+0.001)>=curTh);
                write.table(curRes[,1:4],file=resFile,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE,quote=FALSE);
                curRes_asMatrix <- matrix(unlist(curRes),byrow=FALSE,ncol=5)
#                curSelection <- curRes_asMatrix[which(curRes_asMatrix[,5]=="TRUE"),1:4]
#		positiveWindows <- matrix(unlist(curSelection),byrow=FALSE,ncol=4)
		positiveWindows <- curRes[i,1:4]
                write.table(positiveWindows,file=selectionFile,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE);
		return(positiveWindows)
        })
        print(class(curStatus))
	print(names(x))
	names(curStatus)=names(x)

	print(names(curStatus))
	return(curStatus)
}


#---------------------------
groupPosFrags <- function(data,th)
{
state=0
m=0;n=0;s=0;c="";e=0;t=0

res=c()
for(i in 1:nrow(data))
{
    if(data[i,4]>=th)
    {
        if(state<1){c=as.character(data[i,1]);s=data[i,2];e=data[i,3];n=1;t=data[i,4];m=m+1}
        else{e=data[i,3];n=n+1;t=t+data[i,4]}
        state=1;
    }
    else
    {
        if(state>0){
            res<-rbind(res,cbind(c,s,e,e-s,t/n))
            state=0;c="";s=0;e=0}
    }
}
if(s>0){res<-rbind(res,cbind(c,s,e,e-s,t/n))}
return(res)
}


