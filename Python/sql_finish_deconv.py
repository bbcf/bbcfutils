#!/usr/bin/env python

import rpy2.robjects as robjects
from bbcflib.btrack import Track
import os
import sys

usage = """sql_finish_deconv.py input output

input    input Rdata file
output   output sql file
chrom    chromosome name
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg


def _robject(obj,chrom):
    vals = []
    start = 0
    for p in obj.iter_row():
        if chrom == p.rx2('chr')[0]:
            pos = p.rx2('pos')[0]
            score = p.rx2('score')[0]
            yield((pos-1,pos,score))


def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) != 3:
            raise Usage("sql_finish_deconv.py takes exactly three arguments.")

        infile = argv[0]
        outfile = argv[1]
        chrom = argv[2]
        if not(os.path.exists(infile)):
            raise Usage("File %s does not exist." % infile)
        robjects.r.load('%s' %infile)
        with Track(outfile, format='sql') as track:
            track.write(chrom,_robject(robjects.r.wig,chrom))
        with open(outfile+'_deconv.bed','a') as fbed:
            for p in robjects.r.bed.iter_row():
                bed_row = [p.rx2('chr')[0],
                           str(int(p.rx2('start')[0])),str(int(p.rx2('end')[0])),
                           p.rx2('name')[0],str(p.rx2('score')[0])]
                fbed.write("\t".join(bed_row)+"\n")
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())
