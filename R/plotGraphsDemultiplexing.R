Args <- commandArgs(TRUE)
print(length(Args))
print(Args)
reportFile <- Args[1]
plotFile <- Args[2]

library(RColorBrewer)
library(plotrix)

plotSeparationReport <- function(reportFie)
{
	myColors=c(brewer.pal(9,"Pastel1"),"#99CCFF","#FF9966")

	print(paste("will treat report file",reportFile))
	data_all <- read.delim(reportFile)

	#Daan_run5_Sept2010_lane1_step1.report
	nSeries<-nrow(data_all)-2
	i <- c(1:nSeries)
	data <- data_all[i,]

	#Plots number per primer/unclassified
	par(mfrow=c(2,1))
	data <- data_all[,2]
	names(data) <- data_all[,1]
	dataPercent <- 100*data[1:(nSeries+1)]/data[nrow(data_all)]
	labs <- sprintf("%s\n%.2f%%",names(data)[1:(nSeries+1)],dataPercent)
	pie(as.matrix(data[1:(nSeries+1)]),col=myColors[1:(nSeries+1)],label=labs,main="Separation of reads")
	legend("topright",paste("total number of reads\n",data[nrow(data_all)]),col="white",box.col="white", text.col ="orange")	
	plot(1:10,1:10,col="white",axes=FALSE,xlab="",ylab="")
	data2Plot <- cbind(as.matrix(data[1:(nSeries+1)]))
	colnames(data2Plot) <- c("#reads")
	addtable2plot(x=3.5,y=8,data2Plot,display.rownames=TRUE,display.colnames=TRUE)
	
	#Plots Undigested,... + mappable/Excluded
    	par(mfrow=c(2,1))
  	data <- data_all[i,]
    	dataPercent <- 100*data[,3:4]/data[,2]
   	rownames(dataPercent) <- data[,1]
    	colnames(dataPercent) <- paste("%",colnames(data)[3:4],sep="")

 	barplot(as.matrix(dataPercent[]), col=myColors[1:nSeries],names=c("%Excluded","%mappable"),beside=TRUE,ylab="percent of reads",ylim=c(0,100))
    	legend("topleft",rownames(dataPercent),fill=myColors[1:nSeries])
    	title("#reads Mappable vs. Excluded*")
    	data2Plot <- cbind(as.matrix(data[,3:4]))
    	rownames(data2Plot) <- data[,1]
    	colnames(data2Plot) <- c("Excluded","Mappable")
    	plot(1:10,1:10,col="white",axes=FALSE,xlab="",ylab="")
    	addtable2plot(x=2.5,y=8,data2Plot,display.rownames=TRUE,display.colnames=TRUE)
    	legend("bottomleft",c("* Excluded: reads that have been filtered, due to the presence of Undigested,Self-ligated and/or bait sequence","  Mappable: reads considered for the rest of the analysis (Bowtie)"),cex=0.65)


}


pdf(plotFile)
plotSeparationReport(reportFile)
dev.off()

