#!/archive/epfl/bbcf/bin/bin.x86_64/python
import rpy2.robjects as robjects
import sqlite3
import re
from numpy import *
import getopt
import os
import sys

usage = """sql_prepare_deconv.py sql_fwd sql_rev input_bed output chrname chrlength cutoff read_len

sql_fwd   input sql file for forward strand
sql_rev   input sql file for reverse strand
input_bed input bed file with region to select from
output    output Rdata file
chrname   name of chromosome to fetch from sql input
chrlength length of that chromosome
cutoff    maximum region size to consider
read_len  read length
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def select_bed_line(c,s,e):
    sql = "select max(start,"+str(s)+"), "
    sql += "min(end,"+str(e)+"), "
    sql += "score from '"+c+"' where "
    sql += "end>"+str(s)+" and start<"+str(e)+";"
    return sql

def parse_bed(file_name, seqname=None):
    def select_by_seqname(l):
        if seqname == None:
            return True
        else:
            return re.match(seqname+r'\b',l)
    with open(file_name,'r') as fh:
        bed_list=[l.rstrip('\n').split('\t') for l in fh if select_by_seqname(l)]
        fh.close()
    return bed_list

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) != 8:
            raise Usage("sql_prepare_deconv takes exactly 8 arguments.")

        dbfwd = argv[0]
        dbrev = argv[1]
        bedfile = argv[2]
        output_file = argv[3]
        chromosome_name = argv[4]
        chromosome_length = argv[5]
        size_cutoff = argv[6]
        read_length = argv[7]
#        prod_shift = 50

        strands = {dbfwd:'plus',dbrev:'minus'}
        for db in strands.keys():
            if not(os.path.exists(db)):
                raise Usage("Sqlite file %s does not exist." % db)
        if not(os.path.exists(bedfile)):
            raise Usage("Bed file %s does not exist." % bedfile)
        if os.path.exists(output_file):
            raise Usage("Output file %s already exists." % output_file)
        bed_list = parse_bed(bedfile,seqname=chromosome_name)
        robjects.r('counts=data.frame()')
        row_count = 0
        for bed_row in bed_list:
            row_count += 1
            chr = bed_row[0]
            start = int(bed_row[1])
            end = int(bed_row[2])
            if len(bed_row)>3:
                reg_name = bed_row[3]
            else:
                reg_name = row_count
            if end-start > size_cutoff:
                continue
            if start < 0:
                start = 0
            if end > chromosome_length:
                end = chromosome_length
            allpos = range(start+1,end+1)
            data_block = {'pos': robjects.IntVector(allpos),
                          'plus': robjects.FloatVector([0 for i in allpos]),
                          'minus': robjects.FloatVector([0 for i in allpos]),
                          'name': robjects.StrVector([reg_name for i in allpos])}
            for db,name in strands.iteritems():
                connection = sqlite3.connect(db)
                cur = connection.cursor()
                cur.execute(select_bed_line(chr, start, end))
                connection.commit()
                n = 0
                for sql_row in cur:
                    while data_block['pos'][n] <= sql_row[0]:
                        n+=1
                    for p in range(sql_row[0],sql_row[1]):
                        data_block[name][n] = sql_row[2]
                        n+=1
                cur.close()
            r_block = robjects.DataFrame(data_block)
            robjects.r('counts=rbind(counts,%s)'%r_block.r_repr())
        robjects.r('save(counts,file="%s")' %output_file)
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())

