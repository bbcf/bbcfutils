#!/usr/bin/env python

from bbcflib.common import unique_filename_in
import sys, getopt, os


opts = dict(getopt.getopt(sys.argv[1:],"i:o:n:x:",[])[0])

exportFile=opts['-i']
n=opts.get('-n') or 1
x=opts.get('-x') or 22

print("In fastqToFasta")
print("i="+fqFile)
print("n="+str(n))
print("x="+str(n))

faFile=opts.get('-o') or unique_filename_in()
output=open(faFile,"w")
i=1
n=int(n)
x=int(x)
with open(exportFile,"r") as f:
	for s in f:
        	s=s.strip('\n').split('\t')
                output.write(">line"+str(i)+":"+s[8]+":"+s[9]+":"+s[-1]+s[8][(n-1):(n+x-1)]+"\n")
                i=i+1
output.close()
