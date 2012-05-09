Args <- commandArgs(TRUE)
infile <- Args[1]
curName <- Args[2]
prefName <- Args[3]
regCoord <- Args[4]
wmaxDomainograms <- as.numeric(Args[5])
wmax_BRICKS <- as.numeric(Args[6])
script.path = "/mnt/common/epfl/share"
if (length(Args)>7) {script.path=Args[8]}

source(paste(script.path,"/downstreamAnalysis_Rfunctions.R",sep=''))
source(paste(script.path,"/domainogram_functions.R",sep=''))
#source("/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_SAM/.txt") #for getRegions()
#source("/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_SAM/domainogram_functions.txt")

print(paste("length(Args)=",length(Args)))
print(Args)
if(length(Args)<7){nskip=0}else{nskip=Args[7]}
print(paste("skip=",nskip))

if(is.na(wmaxDomainograms)){wmaxDomainograms=500}
if(is.na(wmax_BRICKS)){wmax_BRICKS=50}


data <- read.delim(as.character(infile),skip=nskip,header=F,stringsAsFactors=F)
if(ncol(data)<4)
{
	data <- read.delim(as.character(infile),skip=nskip,header=F,stringsAsFactors=F,sep="\t")
}


if(length(regCoord)<2)
{

	print("not from profile corrected data")
	print(paste("regCoord=",regCoord))
	regCoordSplit <- unlist(strsplit(gsub("-",":",regCoord,perl=TRUE),":"))
	if(length(regCoordSplit)>=3){
		print(regCoordSplit)
		reg=as.numeric(regCoordSplit[2:3])
		curChr=regCoordSplit[1]
		print(paste("chr to treat=",curChr))
		data <- data[which(data[,1]==curChr),]
		print(paste("Will first split data and call domainograms"))
	        print(regCoordSplit)
	        data_splitted <- getRegions(data,regCoordSplit)
       		resfiles_up=runDomainogram(data_splitted[[2]],paste("domainogram_",curName,"_dataUp",sep=""),wmaxDomainograms,wmax_BRICKS,prefName=paste(prefName,"_dataUp",sep=""))
       		resfiles_down=runDomainogram(data_splitted[[3]],paste("domainogram_",curName,"_dataDown",sep=""),wmaxDomainograms,wmax_BRICKS,prefName=paste(prefName,"_dataDown",sep=""))
		resfiles=c(resfiles_up,resfiles_down)
	}else{
		curChr=regCoord
		print(paste("chr to treat=",curChr))
	        data <- data[which(data[,1]==curChr),]
	        print(paste("prefName"=prefName))
 		resfiles=runDomainogram(data,paste("domainogram_",curName,sep=""),wmaxDomainograms,wmax_BRICKS,prefName=prefName)
	}
	
}else{
	print("from profile corrected data")
	print(paste("prefName"=prefName))
	resfiles=runDomainogram(data,paste("domainogram_",curName,sep=""),wmaxDomainograms,wmax_BRICKS,prefName=prefName)
}

print("####resfiles####")
print(paste(resfiles,collapse="\n"))

write.table(print(paste("####resfiles####:",paste(resfiles,collapse="\n"),sep="\n")),file=paste(prefName,".log",sep=""),quote=FALSE)

# return: list of created files
