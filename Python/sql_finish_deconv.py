#!/archive/epfl/bbcf/bin/bin.x86_64/python
import rpy2.robjects as robjects
import sqlite3
import re
from numpy import *
import getopt
import os

usage = """sql_finish_deconv.py input output chrname

input     input Rdata file
output    output sql file
chrname  chromosome name (sql table to insert into)
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) != 3:
            raise Usage("sql_finish_deconv.py takes exactly three arguments.")

        infile = argv[0]
        outfile = argv[1]
        chrname = argv[2]
        if not(os.path.exists(infile)):
            raise Usage("File %s does not exist." % infile)
        if os.path.exists(outfile):
            raise Usage("File %s already exists." % outfile)
        robjects.r.load('%s' %infile)
        connection = sqlite3.connect( outfile )
        vals = []
        last_score
        start = 0
        for p in robjects.r.wig.iter_row():
            pos = p.rx2('pos')[0]
            vals.append((pos-1,pos,p.rx2('score')[0]))
        connection.executemany('insert into %s (start,end,score) values (?,?,?)'%chrname,vals)
        connection.commit()
        connection.close()
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())

