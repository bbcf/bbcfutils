#llNames.s# --args $fragsFiles combinedFile.txt $names 4,5 $reffile
## input a list of files to be combined, separated by ";".
# files will be grouped by chr, start and end coordinates of fragments (default in columns 1,2,3).
# if a reference file is given, as a bed file with all frags to be considered, then the output will contain all/only those frags
## + a list of names for the columns (if none, then a generic one will be generated)
## an optional idCols list of columns to be merged can be provided (e.g., idCols=4,6 or simply idCols=4). Will be the same for all inputFiles

options(stringsAsFactors = FALSE)
options(scipen=999) # disable the scientific notation

Args <- commandArgs(TRUE)
print("length(Args)=",length(Args))
print("combineTracks.R --args ")
print(Args)
fragsFiles.s <- Args[1]; print(paste("fragsFiles=",fragsFiles.s))
combinedFile <- Args[2]; print(paste("resfile=",combinedFile))
if(length(Args)>2 & nchar(Args[3])>2 ){allNames.s=Args[3]}else{allNames.s=""}
if(length(Args)>3){idCols=as.character(Args[4])}else{idCols="4"}; print(paste("idCols=",idCols))
if(length(Args)>4){defVal=Args[5]}else{defVal="0"}
if(length(Args)>5){reffile=Args[6]}else{reffile=""}

# might be passed as parameters
nskip=1
curSep="\t"

if(length(reffile)>1){
    refFrags=read.delim(reffile,sep=curSep,skip=nskip,header=FALSE)
    print(dim(refFrags))
}
## if no reffile => prepare refFrags with the union of all coords of fragments found

data.l=list()
allFragsCoords=c()
fragsFiles=unlist(strsplit(fragsFiles.s,","))
nFiles=length(fragsFiles)
print(paste(nFiles," to be merged"))
if(nchar(allNames.s)<2){allNames=paste(rep("sample",nFiles),1:nFiles,sep="_")}else{allNames=unlist(strsplit(allNames.s,","))}
idfile=1 #no need to care about the sample name
for(infile in fragsFiles){
    print(paste("read ",infile))
    data <- read.delim(infile,skip=nskip,header=FALSE,sep=curSep)
    curFragsCoords=apply(data,1,function(x){paste(x[1:3],collapse="\t")})
    allIdCols=as.numeric(unlist(strsplit(idCols,",")))
    data.l[[idfile]]=matrix(ncol=length(allIdCols),nrow=nrow(data))
    for(i in 1:length(allIdCols)){
        curIdCol=allIdCols[i]
        if(is.numeric(data[,curIdCol])){val=round(data[,curIdCol],3)}else{val=data[,curIdCol]}
        data.l[[idfile]][,i]=val
    }
#    data.l[[idfile]]=data[,c(idCols)]
    rownames(data.l[[idfile]])=curFragsCoords
    allFragsCoords=c(allFragsCoords,curFragsCoords)
    print(dim(data.l[[idfile]]))
    print(head(data.l[[idfile]]))
    idfile=idfile+1
}

#names(data.l)=allNames
allFragsCoords=unique(sort(allFragsCoords))
length(allFragsCoords)

if(length(reffile)>1){
    allFragsCoords=apply(refFrags,1,function(x){paste(x[1],x[2],x[3],sep="\t")})
}

if(length(grep("NA",defVal))>0){defVal=NA}else{defVal=as.numeric(defVal)}
res <- matrix(defVal,ncol=length(allNames),nrow=length(allFragsCoords))
colnames(res)=allNames
rownames(res)=allFragsCoords

allIdCols=as.numeric(unlist(strsplit(idCols,",")))
ncolumns=length(allIdCols)
s=seq(1,ncol(res),by=ncolumns)
e=s+(ncolumns-1)
for(i in 1:length(data.l)){
res[rownames(data.l[[i]]),s[i]:e[i]]=data.l[[i]]
}

print("object res:")
print(dim(res))
print(head(res))

allCoords=matrix(unlist(lapply(allFragsCoords,function(x){unlist(strsplit(x,"\t"))})),ncol=3,byrow=TRUE)
o <- order(allCoords[,1],as.numeric(allCoords[,2]),as.numeric(allCoords[,3]))
curTable=cbind(allCoords[o,1:3],res[o,])
colnames(curTable)=c("chrFrags","chrStart","chrEnd",allNames)
rownames(curTable)=rownames(res)[o]

print(combinedFile)

write.table(curTable,file=combinedFile,sep="\t",row.names=FALSE,quote=FALSE)

print("Done!!")


