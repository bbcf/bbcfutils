#!/archive/epfl/bbcf/bin/bin.x86_64/python
import rpy2.robjects as robjects
import sqlite3
import getopt
import os
import sys

usage = """sql_finish_deconv.py input output

input    input Rdata file
output   output sql file
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) != 2:
            raise Usage("sql_finish_deconv.py takes exactly two arguments.")

        infile = argv[0]
        outfile = argv[1]
        if not(os.path.exists(infile)):
            raise Usage("File %s does not exist." % infile)
        robjects.r.load('%s' %infile)
        connection = sqlite3.connect( outfile )
        vals = []
        start = 0
        for p in robjects.r.wig.iter_row():
            chr = p.rx2('chr')[0]
            pos = p.rx2('pos')[0]
            score = p.rx2('score')[0]
            vals.append((pos-1,pos,score))
        connection.executemany('insert into %s (start,end,score) values (?,?,?)'%chr,vals)
        connection.commit()
        connection.close()
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

