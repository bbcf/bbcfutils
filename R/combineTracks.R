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
if(length(Args)>2 & nchar(Args[3])>2 ){allNames.s=Args[3]}else{allNames.s=""}; print(paste("allNames.s=",allNames.s,sep=""))
if(length(Args)>3){idCols=as.character(Args[4])}else{idCols="4"}; print(paste("idCols=",idCols));
if(length(Args)>4){defVal=Args[5]}else{defVal="0"}; print(paste("defVal=",defVal,sep=""))
if(length(Args)>5){out_chromosomes=Args[6]}else{out_chromosomes="NA"};print(paste("out_chromosomes=",out_chromosomes,sep=""))
if(length(Args)>6){regToExclude.s=Args[7]}else{regToExclude.s="NA"} #comma separated list of regToExclude, 1 per file or nothing [can be passed as: regFile1,NA,regFil3,NA,NA if no regToExclude defined];
print(paste("regToExclude.s=",regToExclude.s,sep=""))
if(length(Args)>7){reffile=Args[8]}else{reffile=""}; print(paste("reffile=",reffile,sep=""))

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


if(regToExclude.s!="NA"){
    regToExclude.s=gsub("NA","chrNA:0-0",regToExclude.s)
    regToExclude.l=unlist(strsplit(regToExclude.s,","))
    print(paste("length(regToExclude.l)=",length(regToExclude.l),sep=""))
    if(length(regToExclude.l)<nFiles){rm("regToExclude.l"); print("You did not provide enough regToExclude regions")}
}
## will only treat the regToExist regions if the list regToExclude.l exists: exists("regToExclude.l")

if(nchar(allNames.s)<2){allNames=paste(rep("sample",nFiles),1:nFiles,sep="_")}else{allNames=unlist(strsplit(allNames.s,","))}
idfile=1 #no need to care about the sample name
for(infile in fragsFiles){
    print(paste("read ",infile))
    data <- read.delim(infile,skip=nskip,header=FALSE,sep=curSep)
    curFragsCoords=apply(data,1,function(x){paste(c(x[1],as.numeric(x[2]),as.numeric(x[3])),collapse="\t")})
    allIdCols=as.numeric(unlist(strsplit(idCols,",")))
    data.l[[idfile]]=matrix(ncol=length(allIdCols),nrow=nrow(data))
    if(exists("regToExclude.l")){
        chr_curRegToExclude=unlist(strsplit(regToExclude.l[idfile],split=":"))[1]
        coord_curRegToExclude=as.numeric(unlist(strsplit(unlist(strsplit(regToExclude.l[idfile],split=":"))[2],"-")))
    }else{chr_curRegToExclude="NA";coord_curRegToExclude=c(0,0)} ## if no regToExclude region defined, set to NA and 0 => will return no fragments
    names(coord_curRegToExclude)=c("start","end")

    for(i in 1:length(allIdCols)){
        curIdCol=allIdCols[i]
        if(is.numeric(data[,curIdCol])){val=round(data[,curIdCol],3)}else{val=data[,curIdCol]}
        data.l[[idfile]][,i]=val
        I.exclude=which(data[,1]==chr_curRegToExclude & data[,2]<coord_curRegToExclude["end"] & data[,3]>coord_curRegToExclude["start"])
        data.l[[idfile]][I.exclude,i]=NA
    }

#    data.l[[idfile]]=data[,c(idCols)]
    rownames(data.l[[idfile]])=curFragsCoords
    allFragsCoords=c(allFragsCoords,curFragsCoords)
    print(dim(data.l[[idfile]]))
    print(head(data.l[[idfile]]))
    idfile=idfile+1
}
#save.image("../rdata_combineTrack.RData")

#names(data.l)=allNames
allFragsCoords=unique(sort(allFragsCoords))
length(allFragsCoords)

if(length(reffile)>1){
    allFragsCoords=apply(refFrags,1,function(x){paste(x[1],as.numeric(x[2]),as.numeric(x[3]),sep="\t")})
}


if(out_chromosomes!="NA"){
    out_chromosomes.l=unlist(strsplit(out_chromosomes,","))
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
print(paste("nrow=",nrow(curTable)," ,ncol=",ncol(curTable),sep=""))

if(out_chromosomes!="NA"){
    out_chromosomes.l=unlist(strsplit(out_chromosomes,","))
    I_outchr = which(curTable[,1] %in% out_chromosomes.l)
    curTable = curTable[I_outchr,]
    print(paste("filter chromosomes:",out_chromosomes,"\nremains ",nrow(curTable)," rows"),sep="")
}else{print("No filter applied => all chromosomes returns")}

print(combinedFile)

#save.image("../rdata_combineTrack.RData")
write.table(curTable,file=combinedFile,sep="\t",row.names=FALSE,quote=FALSE)

print("Done!!")


