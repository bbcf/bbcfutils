library(plotrix)
########
# Create&plot domainograms
# find, select and plot BRICKS
#######
plotDomainogram <- function(myScores,wmax,nGrad,myCex=1.5)
{
        r <- rank(myScores)
        Rq <- 1-((r-0.5)/length(myScores))
        myCols <- smoothColors("black",(nGrad-2)/2,"purple",(nGrad-2)/2,"red")
	myBreaks=seq(-0.0001,6,length.out=10)

        #if(length(myScores)>1000){myPch='.'}else{myPch=17}
        myPch=17
        for(w in 1:wmax)
        {
                cat(paste(w,";",sep=""))
                if(w==1){
			p <- unlist(lapply(1:length(q),function(i){w=w;pchisq(q[i], df=2*w, ncp=0,lower.tail=FALSE)})) 
			p[p<10^-6]=10^-6
#			cat(paste("length(p)=",length(p),";",sep=""))
			plot(1:length(p),rep(w*(1/wmax),length(p)),ylim=c(0,1),col=myCols[cut(-log10(p),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)
			prev <- q 
                }
                else
                {
			s <- unlist(lapply(w:length(q),function(i){-2*log(q[i])+prev[i-w+1]}))
                        p <- unlist(lapply(1:length(s),function(i){w=w;pchisq(s[i], df=2*w, ncp=0,lower.tail=FALSE)}))
			p[p<10^-6]=10^-6
			#cat(paste("length(p)=",length(p),";y=",(w*(1/wmax)),";x1=",(w/2),";xp=",((w/2)+length(p)),";",sep=""))
                        points(c(1:length(p))+(w/2),rep(w*(1/wmax),length(p)),col=myCols[cut(-log10(p),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)
			prev <- s
                }
        }
#	print(summary(p))
	print("Done!")
#	plot(0,0,ylim=c(0,1),xlim=c(0,wmax))
gradient.rect(0,0.5,1,0.9,col=myCols,gradient="y",border=NA)
        rect(1,0.5,2,0.5)
        rect(1,0.66,2,0.66)
        rect(1,0.78,2,0.78)
        rect(1,0.9,2,0.9)
        text(5,0.5,label=10^-myBreaks[1],pch=3)
        text(5,0.66,label=format(10^-myBreaks[4],scientific=TRUE),pch=3)
        text(5,0.78,label=format(10^-myBreaks[7],scientific=TRUE),pch=3)
        text(5,0.9,label=10^-myBreaks[10],pch=3)
}
#plotDomainogram(dataDown[1:100,4],40,20,1)

plotLegend <- function(imgFile="")
{
	if(nchar(imgFile)==0){X11()}else{print(imgFile);pdf(imgFile)}
	nGrad=10
	myCols <- smoothColors("black",(nGrad-2)/2,"purple",(nGrad-2)/2,"red")		
        plot(0,0,xlim=c(0,5),ylim=c(0.3,1),axes=FALSE,col="white",xlab="",ylab="")
	gradient.rect(0,0.5,1,0.9,col=myCols,gradient="y",border=NA)
        rect(1,0.5,1.2,0.5)
        rect(1,0.66,1.2,0.66)
        rect(1,0.78,1.2,0.78)
        rect(1,0.9,1.2,0.9)
        text(2,0.5,label=1,pch=3)
        text(2,0.66,label=format(10^-2,scientific=TRUE),pch=3)
        text(2,0.78,label=format(10^-4,scientific=TRUE),pch=3)
        text(2,0.9,label=10^-6,pch=3)
	if(nchar(imgFile)>0){dev.off()}
}

# optimized (per line by using the precendent + won't plot the domainogram if
# imgFile=NA
createDomainogram <- function(myScores,wmax=20,nGrad=10,myCex=1.5,imgFile="",decrease=TRUE)
{
        r <- rank(myScores)
	if(decrease){q <- 1-((r-0.5)/length(myScores))}else{q <- (r-0.5)/length(myScores)}
        myCols <- smoothColors("black",(nGrad-2)/2,"purple",(nGrad-2)/2,"red")
        myBreaks=seq(-0.0001,6,length.out=10)
        myPch=17;
        m=matrix(nrow=wmax,ncol=length(myScores))
        print(dim(m))
        if(is.na(imgFile)){print("Will not plot the resulting domainogram")}else{if(nchar(imgFile)>0){print("will plot the resulting domainogram:");print(imgFile);plotLegend(paste(imgFile,"_legend.pdf",sep=""));pdf(imgFile);}else{X11()}}
        for(w in 1:wmax)
        {
                if(w/100-round(w/100)==0){cat(paste(w,";",sep=""))}
                if(w==1){
                        s<-unlist(lapply(1:length(q),function(i){-2*log(q[i])}))
                        p <- unlist(lapply(1:length(s),function(i){w=w;pchisq(s[i],df=2*w,ncp=0,lower.tail=FALSE)}));p[which(p<10^-6)]=10^-6; m[w,]=p; prev <- s
		}
                else
                {
                        s <- unlist(lapply(w:length(q),function(i){-2*log(q[i])+prev[i-w+1]}))
                        p <- unlist(lapply(1:length(s),function(i){w=w;pchisq(s[i], df=2*w, ncp=0,lower.tail=FALSE)}))
                        p[which(p<10^-6)]=10^-6;m[w,w:ncol(m)]=p; prev <- s;
                 }

		if(!is.na(imgFile))
		{
                	if(w==1){plot(c(1:ncol(m)),rep(w*(1/wmax),ncol(m)),xlim=c(0,length(myScores)),ylim=c(0,1),col=myCols[cut(-log10(m[w,]),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)}
                	else{points(c(1:length(m[w,w:ncol(m)]))+(w/2),rep(w*(1/wmax),length(m[w,w:ncol(m)])),col=myCols[cut(-log10(m[w,w:ncol(m)]),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)}
		}
        }
        if(!is.na(imgFile) && nchar(imgFile)>0){dev.off()}
        return(m)
}

createDomainogram_vJ <- function(dataToTreat,wmax=20,nGrad=10,myCex=1.5,imgFile="",decrease=TRUE,main='')
{
    myScores=dataToTreat[,4]
    lensc=length(myScores)
    r <- rank(myScores)
    if(decrease){
        q <- -2*log(1-((r-0.5)/lensc))
    } else {
        q <- -2*log((r-0.5)/lensc)
    }
    myCols <- smoothColors("black",(nGrad-2)/2,"purple",(nGrad-2)/2,"red")
    myBreaks=seq(-0.0001,6,length.out=10)
    myPch=17;
    mflat = q
    prev = q
    for(w in 2:wmax) {
        if(w %% 100==0){cat(paste(w,";",sep=""))}
        prev = q[-(1:(w-1))]+prev[1:(lensc-w+1)]
        mflat=c(mflat,prev)
    }
#    mflat =
#    pchisq(mflat,df=rep(2*1:wmax,lensc+1-1:wmax),ncp=0,lower.tail=FALSE)
    s0 = data.frame(s=round(mflat,digits=6),df=rep(2*1:wmax,lensc+1-1:wmax))
    su = unique(s0)
    mflat = data.frame(pv=pchisq(su$s,df=su$df,ncp=0,lower.tail=FALSE),
      row.names=paste(su$s,su$df,sep='_'))[paste(s0$s,s0$df,sep='_'),1]
    mflat[mflat<1e-6]=1e-6
    wptr = lensc
    print("Make matrice from mflat")
    print(paste("lensc=",lensc))
    print(paste("nrow=wmax=",wmax))
    m = matrix(nrow=wmax,ncol=lensc)
    m[1,] = mflat[1:lensc]
    for (w in 2:wmax)  {
        m[w,w:lensc] = mflat[wptr+1:(lensc+1-w)]
        wptr = wptr+lensc+1-w
    }

     plotDomainogram_vJ(dataToTreat,m,wmax,lensc,myCols,imgFile,plotAxis=TRUE)
#    plotDomainogram_vJ(dataToTreat,m,wmax,lensc,myCols,gsub(".pdf","_noAxis.pdf",imgFile),plotAxis=FALSE)

#	plot(xcoord-rep(1:wmax/2,lensc+1-1:wmax),
#             rep(1:wmax/wmax,lensc+1-1:wmax),
#             ylim=c(0,1),
#             xlab='Position',ylab='Size',
#             main=main, yaxt='n',
#             col=myCols[cut(-log10(mflat),myBreaks,labels=FALSE)],
#             pch=myPch,cex=myCex)
#        axis(side=2,at=seq(0,wmax,length.out=10)/wmax,
#             labels=round(seq(0,wmax,length.out=10),0),las=2)
    return(m)
}

plotDomainogram_vJ <- function(dataToTreat,domainogram.m,wmax,lensc,myCols,imgFile='',plotAxis=TRUE,main='')
{
   m=domainogram.m
   if (is.na(imgFile)) {
        print("Will not plot the resulting domainogram")
    } else {
        if(nchar(imgFile)>0){
            print("Will plot the resulting domainogram in file:")
            print(imgFile)
            plotLegend(paste(imgFile,"_legend.pdf",sep=""))
            pdf(imgFile,width = 7, height = 5.5, pointsize = 8)
            #png(imgFile,width = 1024, height = 800, units = "px")
           # pdf(gsub("pdf","pdf",imgFile),width = round(1024/72), height = round(800/72)) #square of side 1/72 inch, which is also the ‘pixel’ size assumed for graphics
         #   pdf(gsub("pdf","pdf",imgFile),width = 18, height = 13, res=300) #square of side 1/72 inch, which is also the ‘pixel’ size assumed for graphics
	   #  postscript(gsub("pdf","ps",imgFile),width=round(1024/72), height=round(800/72))
	    #jpeg(gsub("pdf","jpeg",imgFile),width=1024, height=800, units = "px", quality=90, bg="black")
        } else {
            X11(width = 12, height = 9)
        }
        xcoord=1:lensc
        for (w in 2:wmax) xcoord=c(xcoord,w:lensc)
        image(1:lensc, 1:wmax/wmax, -log10(t(m)),
	   xlab='position (in Mb)', ylab='size',
           main=main, yaxt='n', xaxt='n', col=myCols )
	if(plotAxis){
        	ticks=seq(0,wmax,length.out=11)[-1]
        	axis(side=2,at=ticks/wmax,labels=as.integer(ticks),las=2)
		I=seq(1,lensc,by=2000)
		ticks=dataToTreat[I,2]+(dataToTreat[I,3]-dataToTreat[I,2])
		ticks=round(ticks/10^6,1)
		axis(side=1,at=I,labels=ticks,las=3)
		#title(sub="position (in Mb)")
	}
	if(nchar(imgFile)>0) dev.off()
    }
}




# optimized (per line by using the precendent
createDomainogram_v1 <- function(myScores,wmax=20,nGrad=10,myCex=1.5,imgFile="")
{
        r <- rank(myScores)
        q <- 1-((r-0.5)/length(myScores))
        myCols <- smoothColors("black",(nGrad-2)/2,"purple",(nGrad-2)/2,"red")
        myBreaks=seq(-0.0001,6,length.out=10)
        myPch=17;
        m=matrix(nrow=wmax,ncol=length(myScores))
        print(dim(m))
        if(nchar(imgFile)>0){print(imgFile);plotLegend(paste(imgFile,"_legend.pdf",sep=""));pdf(imgFile);}else{X11()}
        for(w in 1:wmax)
        {
                if(w/100-round(w/100)==0){cat(paste(w,";",sep=""))}
#                for(i in w:length(myScores))
#                {	
			if(w==1){
			s<-unlist(lapply(1:length(q),function(i){-2*log(q[i])}))
			p <- unlist(lapply(1:length(s),function(i){w=w;pchisq(s[i], df=2*w,ncp=0,lower.tail=FALSE)}));p[which(p<10^-6)]=10^-6; m[w,]=p; prev <- s } 
			#if(w==1){s<--2*log(q[i]); p <- pchisq(s, df=2*w, ncp=0,lower.tail=FALSE); m[w,i]=max(p,10^-6);
                        else
                        {
                               # s <- -2*log(prod(q[(i-w+1):i])); p<-pchisq(s,df=2*w, ncp=0,lower.tail=FALSE); m[w,i]=max(p,10^-6);
				s <- unlist(lapply(w:length(q),function(i){-2*log(q[i])+prev[i-w+1]}))
                        	p <- unlist(lapply(1:length(s),function(i){w=w;pchisq(s[i], df=2*w, ncp=0,lower.tail=FALSE)}))
				p[which(p<10^-6)]=10^-6;m[w,w:ncol(m)]=p; prev <- s;
                                #if(is.na(m[w,i])){print(paste("w=",w,";i=",i,"=>",m[w,i]));print(head(s));print(head(p));return(s);exit}
                        }
 #               }
                if(w==1){plot(c(1:ncol(m)),rep(w*(1/wmax),ncol(m)),xlim=c(0,length(myScores)),ylim=c(0,1),col=myCols[cut(-log10(m[w,]),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)}
                else{points(c(1:length(m[w,w:ncol(m)]))+(w/2),rep(w*(1/wmax),length(m[w,w:ncol(m)])),col=myCols[cut(-log10(m[w,w:ncol(m)]),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)
                #tmp=myCols[cut(-log10(m[w,w:ncol(m)]),myBreaks,labels=FALSE)];if("#FFFFFF" %in% tmp){print(which(tmp=="#FFFFFF"));print("exit cause White %color");return(m)}
                }
        }
        if(nchar(imgFile)>0){dev.off()}
        return(m)
}


createDomainogram_v0 <- function(myScores,wmax=20,nGrad=10,myCex=1.5,imgFile="")
{
        r <- rank(myScores)
        q <- 1-((r-0.5)/length(myScores))
        myCols <- smoothColors("black",(nGrad-2)/2,"purple",(nGrad-2)/2,"red")
        myBreaks=seq(-0.0001,6,length.out=10)
        myPch=17;
	m=matrix(nrow=wmax,ncol=length(myScores)) 
	print(dim(m))
	if(nchar(imgFile)>0){print(imgFile);plotLegend(paste(imgFile,"_legend.pdf",sep=""));pdf(imgFile);}else{X11()}
	for(w in 1:wmax)
	{	
		if(w/100-round(w/100)==0){cat(paste(w,";",sep=""))}
        	for(i in w:length(myScores))
        	{
                	if(w==1){s<--2*log(q[i]); p <- pchisq(s, df=2*w, ncp=0,lower.tail=FALSE); m[w,i]=max(p,10^-6); 
			}
                	else
                	{
				s <- -2*log(prod(q[(i-w+1):i])); p<-pchisq(s, df=2*w, ncp=0,lower.tail=FALSE); m[w,i]=max(p,10^-6); 
				#if(is.na(m[w,i])){print(paste("w=",w,";i=",i,"=>",m[w,i]));print(head(s));print(head(p));return(s);exit}
                	}
        	}
		if(w==1){plot(c(1:ncol(m)),rep(w*(1/wmax),ncol(m)),xlim=c(0,length(myScores)),ylim=c(0,1),col=myCols[cut(-log10(m[w,]),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)}
		else{points(c(1:length(m[w,w:ncol(m)]))+(w/2),rep(w*(1/wmax),length(m[w,w:ncol(m)])),col=myCols[cut(-log10(m[w,w:ncol(m)]),myBreaks,labels=FALSE)],pch=myPch,cex=myCex)
		#tmp=myCols[cut(-log10(m[w,w:ncol(m)]),myBreaks,labels=FALSE)];if("#FFFFFF" %in% tmp){print(which(tmp=="#FFFFFF"));print("exit cause White color");return(m)}
		}
	}
	if(nchar(imgFile)>0){dev.off()}
	return(m)
}
#res<-createDomainogram(myScores[1:500],500,10,0.7)

createBRICK <-function(P,wmax,gamma=0.01)
{
#	min_p<-which(P==min(P,na.rm=TRUE))
#	length_min_p=which(P==min_p/nrow(P)
	K=nrow(P);I=ncol(P)
	min_p=0;length_min_p=NA;V=rep(NA,I);Vl<-rep(NA,I);Id=rep(NA,I)
	if(wmax>nrow(P)){print("Error, wmax>nrow(P)");return(NA)}
	for(i in 1:I)
	{	w=1;min_p=0;length_min_p=NA;indice_min_p=NA
		while(w<= min(i,wmax))
		{
			if(w==1){np=gamma*log(P[w,i])}else{np=log(P[w,i])}
			np=w*np; if(i-w<=0){curV=0}else{curV=V[i-w]}; np=np+curV;
			if(np==0){cat(paste(i,",",w,"=>",P[w,i],";",sep="")) }
			if(np<min_p){min_p=np;length_min_p=w;indice_min_p=i}
			w=w+1
		}
		V[i]=min_p;Vl[i]=length_min_p;Id[i]=indice_min_p;
		if(i<4){print("next i")}
	}
	print(paste("Last i=",i,";last w=",w,"=>V[i]=",V[i]))
	print("Traceback"); 
	w=I;n=0;b<-c();
	while(w>0){l=Vl[w];
		b<-rbind(b,c(l,V[w],Id[w],P[l,Id[w]]))
		#print(paste("Vl[",w,"]=>","length=",l,";p=",V[w],";",Id[w],";",sep=""));
		w=w-l;n=n+1}
	print(paste("#BRICKS found=",n))
	return(b)
}
#foundBricks<-createBRICK(res,50)

getBRICKSCoord <- function(foundBRICKS,data,domainogram)
{
	i_start<-foundBRICKS[,3]-foundBRICKS[,1]+1
	i_end<-foundBRICKS[,3]
	avgScores<-as.vector(unlist(lapply(1:length(i_start),function(j){mean(data[i_start[j]:i_end[j],4])})))
	coordBRICKS<-cbind(data[i_start,2],data[i_end,3],foundBRICKS[,1],foundBRICKS[,4],foundBRICKS[,3],avgScores)
	pval <- unlist(lapply(1:nrow(coordBRICKS),function(i){domainogram[coordBRICKS[i,3],coordBRICKS[i,5]]}))
	coordBRICKS<-cbind(coordBRICKS,pval,coordBRICKS[,2]-coordBRICKS[,1])
	colnames(coordBRICKS) <- c("start","end","w_length","BRICKS_pval","idLastFrag","avgFragsScores","Segment_pval","Segment_length")
	o<-order(coordBRICKS[,1]); coordBRICKS<-coordBRICKS[o,]
	return(coordBRICKS)
}
#coordBricks<-getBRICKSCoord(foundBricks,dataDown[1:100,],res)

plotSelectedBRICKS <- function(coordBRICKS,data,myLim=0,myTitle="",imgFile="")
{
	print(paste("******plotSelectedBRICKS*******"))
	i <- which(abs(coordBRICKS[,6])>0)
	j <- which(coordBRICKS[,3]>0)
	k<-intersect(j,intersect(i,which(-log10(coordBRICKS[,4])>=2)))
	print(paste("length(i)=",length(i),"length(j)=",length(j),"length(k)=",length(k),sep="")); print(table(coordBRICKS[k,3]))
	if(length(k)==0){k=1:nrow(coordBRICKS)}
	if(nchar(imgFile)>0){print(imgFile);pdf(imgFile,width=4, height=7)}
	#if(nchar(imgFile)>0){print(imgFile);png(imgFile,width=1024, height=1800)}
	myCols <- smoothColors("black",max(coordBRICKS[,3])-2,"red")
	#myBreaks=seq(-0.0001,6,length.out=10)

	#colToUse=myCols[cut(coordBRICKS[k,3],myBreaks,labels=FALSE)]
	colToUse=rainbow(length(names(table(coordBRICKS[k,3]))))

	par(mfrow=c(3,1))
	if(myLim==0){
		plot(data[,c(2,4)],xlim=c(min(data[,2]),max(data[,2])),ylim=c(1.2*min(data[,4]),1.2*max(data[,4])),type="h",col="grey")
	}
	else{
		plot(data[,c(2,4)],xlim=c(min(data[,2]),min(data[,2])+myLim),ylim=c(1.2*min(data[,4]),1.2*max(data[,4])),type="h",col="grey")
	}
	rect(coordBRICKS[k,1],max(data[,4]),coordBRICKS[k,2],1.2*max(data[,4])+coordBRICKS[k,3],col=colToUse,border=colToUse)
	abline(h=0)

        if(myLim==0){
		plot(c(0,0),col="white",xlim=c(min(coordBRICKS[,1]),max(coordBRICKS[,2])),ylim=c(1.2*min(-log10(coordBRICKS[k,4])),1.2*max(-log10(coordBRICKS[,4]))),xlab="",ylab="-log10(p-val)",main="pval of selected BRICKS")
	}
	else{
		plot(c(0,0),col="white",xlim=c(min(coordBRICKS[,1]),min(coordBRICKS[,2])+myLim),ylim=c(0,1.2*max(-log10(coordBRICKS[,4]))),xlab="",ylab="-log10(p-val)",main="pval of selected BRICKS")

	}
	rect(coordBRICKS[k,1],-log10(coordBRICKS[k,4]),coordBRICKS[k,2],-log10(coordBRICKS[k,4])+0.2,col=colToUse,border=colToUse)
	
#	myBy=5
#	grps=cut(coordBRICKS[k,3],breaks=seq(0,max(coordBRICKS[k,3])+1,by=myBy),right=TRUE)
#	ngrps=length(seq(0,max(coordBRICKS[k,3])+1,by=myBy))
	if(myLim==0){
                plot(c(0,0),col="white",xlim=c(min(coordBRICKS[,1]),max(coordBRICKS[,2])),ylim=c(min(coordBRICKS[k,3])-1,max(coordBRICKS[,3])),xlab="",ylab="segment length",main="length of selected BRICKS",axes=TRUE)
#		valX=seq(0,round(max(coordBRICKS[,1])/10^7),by=2)
#		plot(c(0,0),col="white",xlim=c(min(coordBRICKS[,1]),max(coordBRICKS[,2])),ylim=c(0.4,0.5*ngrps),xlab="",ylab="segment length",main="length of selected BRICKS",axes=FALSE)
		#axis(1,at=10^7*valX,labels=valX)
#		axis(2,at=0.1+0.5*1:(ngrps-1),labels=levels(grps),las=1,cex=0.8)
		#axis(2,at=1:length(names(table(coordBRICKS[,3]))),labels=names(table(coordBRICKS[,3])))
        }
        else{
                plot(c(0,0),col="white",xlim=c(min(coordBRICKS[,1]),min(coordBRICKS[,2])+myLim),ylim=c(0,max(coordBRICKS[,3])+1),xlab="",ylab="segment length",main="length of selected BRICKS")

        }
#        rect(coordBRICKS[k,1],coordBRICKS[k,3],coordBRICKS[k,2],coordBRICKS[k,3]+0.2,col=colToUse,border=colToUse)
        rect(coordBRICKS[k,1],coordBRICKS[k,3],coordBRICKS[k,2],coordBRICKS[k,3]+0.2,col=colToUse,border=colToUse)
#	rect(coordBRICKS[k,1],0.5*as.numeric(grps),coordBRICKS[k,2],0.5*as.numeric(grps)+0.2,col=colToUse,border=colToUse)

	if(nchar(imgFile)>0){dev.off()}
	if(length(k)==1){print("k=1=>will force the res to be a matrix");return(matrix(coordBRICKS[k,],nrow=k,ncol=ncol(coordBRICKS),byrow=TRUE))}
	else{return(coordBRICKS[k,])}
}

getAnnotationsBRICKS <-function(selectedBRICKS,extension=0,curName="",specie="mouse",curChr="2",path=NA)
{
library(ChIPpeakAnno)

#if(specie=="human"){data(TSS.human.GRCh37);myTSS=TSS.human.GRCh37}else{data(TSS.mouse.NCBIM37);myTSS=TSS.mouse.NCBIM37}
if(specie=="human"){
	ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
}
else
{
	ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
}
options(stringsAsFactors=F)
library(biomaRt)

if(is.na(path)){path=""}
# !! chromosome name defined here (in space)
regions = RangedData(IRanges(start=selectedBRICKS[,1]-extension,end=selectedBRICKS[,2]+extension,names=paste("selectedBRICKS",1:nrow(selectedBRICKS),sep="_")),space=rep(curChr,nrow(selectedBRICKS)))
#annotatedBricks <- annotatePeakInBatch(regions, AnnotationData=myTSS)
annotatedBricks <- annotatePeakInBatch(regions, mart=ensembl, featureType=c("TSS"))
annotatedBricks.df <- as.data.frame(annotatedBricks)
i_haveOverlap <- which(annotatedBricks.df$insideFeature %in% c("includeFeature","inside","overlapEnd","overlapStart"))
o <- order(as.numeric(annotatedBricks.df$start))
write.table(annotatedBricks.df[o,],paste(path,curName,"_selectedBRICKS_withAnno.txt",sep=""),row.names=FALSE,quote=FALSE,sep="\t")

#library(topGO)
options(stringsAsFactors=F)
#ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
#att.sel=c("refseq_dna","ensembl_gene_id","mgi_symbol","chromosome_name","start_position","end_position","description")
foundGenes <- getGene( annotatedBricks.df$feature, type="ensembl_gene_id", ensembl)

header=paste("track type=bed name='",curName," (genes in BRICKS)' description='",curName," (genes in BRICKS - extension of +/-",round(extension/1000),"kb )' visibility=dense",sep="")
write.table(header,paste(path,curName,"_foundGenes.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
strands <- foundGenes$strand; strands[strands<0]="-"; strands[strands>=0]="+"
toWrite <- cbind(paste("chr",foundGenes$chromosome_name,sep=""),foundGenes$start_position,foundGenes$end_position,paste(foundGenes$ensembl_gene_id,foundGenes$mgi_symbol,sep=";"),rep(1,nrow(foundGenes)),strands)
write.table(toWrite,paste(path,curName,"_foundGenes.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}


runDomainogram <- function(dataToTreat,curName,wmaxDomainogram=500,wmax_BRICKS=50,decrease=TRUE,specie="mouse",path=NA,myCex=0.7,prefName=NA)
{
        #resDomainogram<-createDomainogram(dataToTreat[,4],wmaxDomainogram,10,0.7,imgFile=paste(curName,".pdf",sep=""))
#        resDomainogram<-createDomainogram(dataToTreat[,4],wmaxDomainogram,10,0.7,imgFile=NA)
	o <- order(dataToTreat[,2])
	dataToTreat=dataToTreat[o,]
	#dataToTreat[dataToTreat[,4]<1,4]=0
	
	curChr=as.character(dataToTreat[1,1])
	
	if(is.na(prefName)){prefName=curName}
	if(is.na(path)){path=""}

	print("In runDomainograms:")
	print(paste("wmaxDomainogram=",wmaxDomainogram,",wmax_BRICKS=",wmax_BRICKS,",decrease=",decrease,",prefName=",prefName))        
	print(head(dataToTreat))
	
	resDomainogram<-createDomainogram_vJ(dataToTreat,wmaxDomainogram,10,myCex=myCex,imgFile=paste(path,prefName,"_domainogram.pdf",sep=""),decrease=decrease)
#        resDomainogram<-createDomainogram_vJ(dataToTreat[,4],wmaxDomainogram,10,0.7,imgFile=NA,decrease=decrease)
	print("createBRICK")
        foundBricks<-createBRICK(resDomainogram,wmax_BRICKS)
        coordBricks<-getBRICKSCoord(foundBricks,dataToTreat,resDomainogram)
        write.table(cbind(rep(curChr,nrow(coordBricks)),coordBricks),paste(path,prefName,"_foundBRICKS.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	save(dataToTreat,resDomainogram,foundBricks,coordBricks,file=paste(path,prefName,"_domainograms.RData",sep=""))
	selectedBricks<-plotSelectedBRICKS(coordBricks,dataToTreat,myLim=0,myTitle=curName,imgFile=paste(path,prefName,"_selectedBRICKS.pdf",sep=""))
        write.table(cbind(rep(curChr,nrow(selectedBricks)),selectedBricks),paste(path,prefName,"_selectedBRICKS.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

	#Write bedGraph selected BRICKS	
	selectedBricksFile=paste(path,prefName,"_selectedBRICKS.bedGraph",sep="")
	header=paste("track type=bedGraph name='",curName," (selected BRICKS)' description='",curName," (selected BRICKS)' visibility=full windowingFunction=maximum",sep="")
	write.table(header,selectedBricksFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
	toWrite<-cbind(rep(curChr,nrow(selectedBricks)),selectedBricks[,1:2],round(selectedBricks[,6],2))
	write.table(toWrite,selectedBricksFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)

	selectedBricksFile_pval=paste(path,prefName,"_selectedBRICKS_pval.bedGraph",sep="")
	header=paste("track type=bedGraph name='",curName," (-log10(pval) selected BRICKS)' description='",curName," (-log10(pval) selected BRICKS)' visibility=full windowingFunction=maximum",sep="")
	write.table(header,selectedBricksFile_pval,row.names=FALSE,col.names=FALSE,quote=FALSE)
	toWrite<-cbind(rep(curChr,nrow(selectedBricks)),selectedBricks[,1:2],round(selectedBricks[,4],2))
	write.table(toWrite,selectedBricksFile_pval,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)

	gsub("chr","",curChr)
	save(dataToTreat,resDomainogram,selectedBricks,foundBricks,coordBricks,file=paste(path,prefName,"_domainograms.RData",sep=""))
# doesn't work when called via lsf ...
#	print(paste("---- Will call getAnnotationsBRICKS with selectedBricks (",nrow(selectedBricks)," BRICKS)",sep=""))
#	if(specie=="human") #pb with "mouse"...
#	{
#		getAnnotationsBRICKS(selectedBricks,extension=1000,curName=curName,specie=specie,curChr=curChr,path=path)
#	}
	print("Done!")
	resFiles=c(paste(path,prefName,"_domainogram.pdf",sep=""),paste(path,prefName,"_foundBRICKS.txt",sep=""),paste(path,prefName,"_domainograms.RData",sep=""),paste(path,prefName,"_selectedBRICKS.pdf",sep=""),paste(path,prefName,"_selectedBRICKS.txt",sep=""),selectedBricksFile,selectedBricksFile_pval)

	return(resFiles)
}


