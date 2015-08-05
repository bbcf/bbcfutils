#mergeRep.R --args infiles mergefile idColScore regToExclude
# infiles="file1;file2;..."  with coordinates are in columns 1,2,3. For now the separator is "\t" and skip=1
# if at least 3 parameters, scores are in column idColScore (default=4)
# if at least 4 parameters => a region to be excluded is provided, otherwise it is set to a fake one which does not affect the process.
# e.g.,
# mergeRep.R --args file1;file2;file3 mergedfile.bedGraph 4 chr2:73911294-73916104

options(stringsAsFactors = FALSE)
options(scipen=999) # disable the scientific notation

Args <- commandArgs(TRUE)
print(length(Args))
print("mergeRep.R --args ")
print(Args)
stopifnot(length(Args)>2)

fragsFiles.s <- Args[1]
mergeFile <- Args[2]
if(length(Args)>2){idColScore=as.numeric(Args[3])}else{idColScore=4}
if(length(Args)>3 & nchar(Args[4])>2){regToExclude=Args[4]}else{regToExclude="chr2:1000-1000"} #default is a fake region which does not concern any fragments

print(paste("fragsFiles.s=",fragsFiles.s,sep=""))
print(paste("mergeFile=",mergeFile,sep=""))
print(paste("idColScore=",idColScore,sep=""))
print(paste("regToExclude=",regToExclude,sep=""))


# might be passed as parameters
nskip=1
curSep="\t"

chr_regToExclude=unlist(strsplit(regToExclude,split=":"))[1]
coord_regToExclude=as.numeric(unlist(strsplit(unlist(strsplit(regToExclude,split=":"))[2],"-")))
names(coord_regToExclude)=c("start","end")


data.l=list()
allFragsCoords=c()
fragsFiles=unlist(strsplit(fragsFiles.s,","))
idfile=1 #no need to care about the sample name
for(infile in fragsFiles){
    print(infile)
    data <- read.delim(infile,skip=nskip,header=FALSE,sep=curSep)
    colnames(data)=c("chromosome","start","end","score")
    curFragsCoords=apply(data,1,function(x){paste(c(x[1],as.numeric(x[2]),as.numeric(x[3])),collapse="\t")})
    data.l[[idfile]]=as.matrix(data[,idColScore])
    rownames(data.l[[idfile]])=curFragsCoords
    colnames(data.l[[idfile]])=infile
    allFragsCoords=c(allFragsCoords,curFragsCoords)
    I.exclude=which(data[,"chromosome"]==chr_regToExclude & data[,"start"]<coord_regToExclude["end"] & data[,"end"]>coord_regToExclude["start"])
    print(length(I.exclude));
    data.l[[idfile]][I.exclude,1]=NA
    idfile=idfile+1
}

names(data.l)=1:length(data.l)
allFragsCoords=unique(sort(allFragsCoords))
length(allFragsCoords)

res <- matrix(0,ncol=length(data.l),nrow=length(allFragsCoords))
colnames(res)=names(data.l)
rownames(res)=allFragsCoords

for(curName in names(data.l)){
res[rownames(data.l[[curName]]),curName]=data.l[[curName]][,1]
}


res.mean=rowMeans(res,na.rm=TRUE)
allCoords=matrix(unlist(lapply(names(res.mean),function(x){unlist(strsplit(x,"\t"))})),ncol=3,byrow=TRUE)
res.mean.df=data.frame(chromosome=allCoords[,1],start=as.numeric(allCoords[,2]), end=as.numeric(allCoords[,3]), mean.score=round(res.mean,2))
o <- order(res.mean.df[,1], res.mean.df[,2],res.mean.df[,3])


write.table(res.mean.df[o,],file=mergeFile,sep="\t",quote=FALSE,row.names=FALSE)

print("Done!!")

