Args <- commandArgs(TRUE)
print(length(Args))
print(Args)
reportFile <- Args[1]
plotFile <- Args[2]

options(stringsAsFactors = FALSE)
library(RColorBrewer)
library(plotrix)

plotSeparationReport <- function(reportFie)
{
    myColors=c(brewer.pal(9,"Pastel1"),"#99CCFF","#FF9966")

    print(paste("will treat report file",reportFile))
    data_all <- read.delim(reportFile)

    #Daan_run5_Sept2010_lane1_step1.report
    I=which(!data_all[,1] %in% c("Total","Unclassified","Discarded","Ambiguous"))
    nSeries<-length(I)
    i <- c(1:nSeries)
    data <- data_all[i,]
    i_ambiguous=which(data_all[,1]=="Ambiguous")
    i_discarded=which(data_all[,1]=="Discarded")
    i_unclassified=which(data_all[,1]=="Unclassified")
    i_tot=which(data_all[,1]=="Total")

    tot_unclassified=sum(data_all[c(i_ambiguous,i_discarded,i_unclassified),2])
    tot=data_all[i_tot,2]

    #Plots number per primer/unclassified
    par(mfrow=c(2,1))
    data <- data_all[,2]
    names(data) <- data_all[,1]
#    dataPercent <- 100*data[1:(nSeries+1)]/data[nrow(data_all)]
    dataPercent <- 100*data/data[nrow(data_all)]
    dataToPlot <- as.matrix(cbind(c(data[I],tot_unclassified),c(dataPercent[I],100*tot_unclassified/tot)))
    rownames(dataToPlot)[nrow(dataToPlot)]="Unclassified*"
    labs <- sprintf("%s\n%.1f%%",rownames(dataToPlot),dataToPlot[,2])
    pie(dataToPlot,col=myColors,label=labs,main="Separation of reads",cex=0.8)
    legend("topright",paste("total number of reads\n",formatC(tot,big.mark = ",", format = "d")),col="white",box.col="white", text.col ="orange")
    plot(1:10,1:10,col="white",axes=F,xlab="",ylab="");
    o <- order(dataToPlot[,1])
    data2Plot=cbind(formatC(dataToPlot[o,1],big.mark = ",", format = "d"),paste(format(dataToPlot[o,2],digits=2),"%",sep=""));
    colnames(data2Plot) <- c("#reads","%reads")
    addtable2plot(x=3.5,y=5,data2Plot,display.rownames=TRUE,display.colnames=TRUE, hlines = TRUE,  vlines = TRUE,bty="o")
    data2Plot <- cbind(formatC(data_all[c(i_ambiguous,i_discarded),2],big.mark = ",", format = "d"),sprintf("%.1f%%",dataPercent[c("Ambiguous","Discarded")]))
    rownames(data2Plot)=c("Ambiguous","Discarded")
    colnames(data2Plot) <- c("#reads","%reads")
    addtable2plot(x=3.5,y=1,data2Plot,display.rownames=TRUE,display.colnames=T, hlines = TRUE,  vlines = TRUE,bty="o",title="*Among Unclassified are:",cex=0.85)
    mtext("Ambiguous: sequences with equally valid classifications\nDiscarded: sequences that are too short (<l) after barcode trimming",side=1,cex=0.6,font=3)

    #Plots Undigested,... + mappable/Excluded
    par(mfrow=c(2,1))
    data <- data_all[i,]
    dataPercent <- 100*data[,3:4]/data[,2]
    rownames(dataPercent) <- data[,1]
    colnames(dataPercent) <- paste("%",colnames(data)[3:4],sep="")
    barplot(as.matrix(dataPercent[]), col=myColors[1:nSeries],names=c("%Excluded","%mappable"),beside=TRUE,ylab="percent of reads",ylim=c(0,100))
    legend("topleft",rownames(dataPercent),fill=myColors[1:nSeries],cex=0.85)
    title("#reads Mappable vs. Excluded*")
    o <- order(data[,4])
    data2Plot <- formatC(cbind(as.matrix(data[o,3:4])),big.mark = ",", format = "d")
    rownames(data2Plot) <- data[o,1]
    colnames(data2Plot) <- c("Excluded","Mappable")
    plot(1:10,1:10,col="white",axes=FALSE,xlab="",ylab="")
    addtable2plot(x=2.5,y=8,data2Plot,display.rownames=TRUE,display.colnames=TRUE, hlines = TRUE,  vlines = TRUE, bty="o")
    legend("bottomleft",c("* Excluded: reads that have been filtered, due to the presence of Undigested,Self-ligated and/or bait sequence","  Mappable: reads considered for the rest of the analysis (Bowtie)"),cex=0.65)

}


pdf(plotFile)
plotSeparationReport(reportFile)
dev.off()

