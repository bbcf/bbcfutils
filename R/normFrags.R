## need baitCoords + regToExclude. Take the mean scores of fragments falling into the baitCoords +/- extSize around the bait, fragments into regToExclude not taken into account.
## return the normalised scores.

options(stringsAsFactors = FALSE)
options(scipen=999) # disable the scientific notation

Args <- commandArgs(TRUE)
print(length(Args))
print("normFrags.R --args")
print(Args)
fragsFile <- Args[1]
normFile <- Args[2]
baitCoords <- Args[3]
extSize <- as.numeric(Args[4])
regToExclude="chr2:1000-1000";
if(length(Args)>4){regToExclude=Args[5]; print(paste("regToExclude=",regToExclude))} #default is a fake region which does not concern any fragments
print(paste("regToExclude=",regToExclude))

# might be passed as parameters
nskip=1
curSep="\t"

chr_regToInclude=unlist(strsplit(baitCoords,split=":"))[1]
coord_baitCoords=as.numeric(unlist(strsplit(unlist(strsplit(baitCoords,split=":"))[2],"-")))
mid_baitCoords=round(sum(coord_baitCoords)/2)
if(mid_baitCoords-extSize>0){regStart=mid_baitCoords-extSize}else{regStart=coord_baitCoords[1]} ## make sure the extended region has a positive start coord
coord_regToInclude=c(regStart,mid_baitCoords+extSize)
names(coord_regToInclude)=c("start","end")

chr_regToExclude=unlist(strsplit(regToExclude,split=":"))[1]
coord_regToExclude=as.numeric(unlist(strsplit(unlist(strsplit(regToExclude,split=":"))[2],"-")))
names(coord_regToExclude)=c("start","end")

data <- read.delim(fragsFile,skip=nskip,header=FALSE,sep=curSep)
colnames(data)=c("chr","start","end","score")
I.include=which(data[,"chr"]==chr_regToInclude & data[,"start"]<coord_regToInclude["end"] & data[,"end"]>coord_regToInclude["start"])
I.exclude=which(data[,"chr"]==chr_regToExclude & data[,"start"]<coord_regToExclude["end"] & data[,"end"]>coord_regToExclude["start"])
data[I.exclude,4]=NA

data.reg=data[I.include,]
mean.reg=mean(data.reg[,4],na.rm=TRUE)
print(mean.reg)


data.norm=data.frame(chr=data[,1], start=as.numeric(data[,2]), end=as.numeric(data[,3]), norm.score=as.numeric(data[,4]/mean.reg))
o <- order(data.norm[,1], data.norm[,2],data.norm[,3])
write.table(data.norm[o,],file=normFile,sep="\t",quote=FALSE,row.names=FALSE)

save.image("../res_normFrags.RData")

print("Done!")

