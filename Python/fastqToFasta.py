#!/usr/bin/env python

from bbcflib.common import unique_filename_in
import sys, getopt, os


def getReadLength(fqFile):
        i=1
        read_length=0
        with open(fqFile,"r") as f:
                for s in f:
                        s=s.strip('\n')
                        if i==2:
                                read_length=len(s)
                        if i>2:
                                return read_length
                        i=i+1
        return read_length



opts = dict(getopt.getopt(sys.argv[1:],"i:o:n:x:",[])[0])

fqFile=opts['-i']

n=opts.get('-n') or 1
x=opts.get('-x') or 22

print("In fastqToFasta")
print("i="+fqFile)
print("n="+str(n))
print("x="+str(n))

faFile=opts.get('-o') or unique_filename_in()
output=open(faFile,"w")
i=1; nextIsQual=0; nextIsSeq=0;
n=int(n);x=int(x)
read_length=getReadLength(fqFile)
print("readLength="+str(read_length))
with open(fqFile,"r") as f:
	for s in f:
        	s=s.strip('\n')
                i=i+1
                if re.search(r'^@',s) and len(s)<read_length:
                	nextIsSeq=1
                if re.search(r'^\+',s) and len(s)<read_length:
                	nextIsQual=1
                        nextIsSeq=0
                if len(s)>=read_length and nextIsSeq>0:
                	seq="".join(s.split('\t')[0][(n-1):(n+x-1)])
                        allSeq=s
                        nextIsSeq=0
                if len(s)>=read_length and nextIsQual>0:
                	qual=s
                        output.write(">line"+str(i)+"_"+allSeq+"_"+qual+"\n"+seq+"\n")
                        nextIsQual=0
#	output.write(">line"+str(i)+"_"+allSeq+"_"+qual+"\n"+seq+"\n")

#return faFile

