BEGIN{th=0.35;tot=0;normFactor=1000;totScore=0;}
{
FS="\t"
n=split($4,a,"|");
if(NR==1)
{
#	print "regToExclude (before split)="reg2Excl;
	nReg=split(reg2Excl,allRegToExcl,",");
#	print "nReg="nReg
	for(idReg in allRegToExcl)
	{
#		print "curRegToExclude(before)="allRegToExcl[idReg];
		gsub(":","\t",allRegToExcl[idReg]);
		gsub("-","\t",allRegToExcl[idReg]);
#		print "curRegToExclude="allRegToExcl[idReg];
	}
}
isToExclude=0
for(idReg in allRegToExcl)
{
	#print "will split coord for reg:"allRegToExcl[idReg]
	split(allRegToExcl[idReg],reg2Excl_coord,"\t");
	#print "split coord=>"reg2Excl_coord[1]";"reg2Excl_coord[2]";"reg2Excl_coord[3];
	if(($1 == reg2Excl_coord[1] && $2 > reg2Excl_coord[2] && $2 < reg2Excl_coord[3]) || ($1 == reg2Excl_coord[1] && $3 > reg2Excl_coord[2] && $3 < reg2Excl_coord[3]))
	{isToExclude=1}
}
if(isToExclude>0)
{nExcluded++;
#print "fragExcluded:"$0;
}
else
{
tot=tot+$5;
#print "$5="$5"=>tot="tot
}
if(frag[$1":"a[2]]){if(a[1]~/endSegment/){endSeg[$1":"a[2]]=a[n]"\t"$5;}else{startSeg[$1":"a[2]]=a[n]"\t"$5;};frag[$1":"a[2]]=frag[$1":"a[2]]"\t"startSeg[$1":"a[2]]"\t"endSeg[$1":"a[2]]}
else{frag[$1":"a[2]]=a[3]"\t"a[5]; if(a[1]~/endSegment/){endSeg[$1":"a[2]]=a[n]"\t"$5;}else{startSeg[$1":"a[2]]=a[n]"\t"$5;};nFrags++;} #!! $5 => score in bed format
}

END{
nFrags=0;
#print "nExcluded="nExcluded;
for(i in frag)
{
nFrags++;
n=split(frag[i],b,"\t");
if(n==4){frag[i]=frag[i]"\t"b[3]"\t"b[4];b[5]=b[3];b[6]=b[4]}
split(b[1],c,":"); split(c[2],d,"-");
if((c[1] == reg2Excl_coord[1] && d[1] > reg2Excl_coord[2] && d[1] < reg2Excl_coord[3]) || (c[1] == reg2Excl_coord[1] && d[2] > reg2Excl_coord[2] && d[2] < reg2Excl_coord[3]))
{
score="NA"; source="ExcludedFrag" #was "_notValid_Excluded"
}
else
{
if(frag[i] ~ /_badStartFrag/)
{
	if(b[5]<th){score=2*b[6]; totScore=totScore+score; source="badStartFrag_but_notRepeatEndFrag";}
	else{score="NA"; source="badStartFrag_and_repeatEndFrag"}
}
else
{
        if(frag[i] ~ /_badEndFrag/)
        {
		if(b[3]<th){score=2*b[4]; totScore=totScore+score; source="badFragEnd_but_notRepeatStartFrag"}
		else{score="NA";source="badEndFrag_and_repeatStartFrag"}
        }
        else
        {
                if(frag[i] ~ /FragIsValid/)
                {
                        if(b[3]<th && b[5]<th){score=b[4]+b[6]; totScore=totScore+score; source="bothValidAndNotRepeated";}
			else{
				if(b[3]<th){score=2*b[4];totScore=totScore+score; source="bothValidButRepeatEndFrag";}
				else{
					if(b[5]<th){score=2*b[6]; totScore=totScore+score; source="bothValidButRepeatStartFrag";}
					else{score="NA";source="bothValidButBothRepeats"}
				}
			}
                }
                else{score="NA"; source="_notValid"}
        }
}
} #end if in regToExclude
print i"\t"frag[i]"\t||\t"source"\t"score"\t"(normFactor*score/tot)"\t"(1000*normFactor*score/tot);
##print i"\t"frag[i]"\t||\t"source"\t"score"\t"(normFactor*score/totScore)"\t"(1000*normFactor*score/totScore);

} #end For
} # end END
