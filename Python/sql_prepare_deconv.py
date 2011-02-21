#!/archive/epfl/bbcf/bin/bin.x86_64/python
import rpy2.robjects as robjects
import sqlite3
import re
from numpy import *
import getopt
import os
import sys

usage = """sql_prepare_deconv.py input_sql input_bed output chrname chrlength cutoff read_len

input_sql input sql files base name (will read input_fwd.sql and inpit_rev.sql)
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
    sql += "score from "+c+" where "
    sql += "end>="+str(s)+" and start<"+str(e)+";"
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
        if len(argv) != 7:
            raise Usage("sql_prepare_deconv takes exactly 7 arguments.")

        db = argv[0]
        bedfile = argv[1]
        output_file = argv[2]
        chromosome_name = argv[3]
        chromosome_length = argv[4]
        size_cutoff = argv[5]
        read_length = argv[6]
#        prod_shift = 50

        strands = {'_fwd.sql':'plus','_rev.sql':'minus'}
        for suffix in strands.keys():
            if not(os.path.exists(db+suffix)):
                raise Usage("Sqlite file %s does not exist." % db+suffix)
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
            if len(bed_row)>2:
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
#                          'prod': robjects.FloatVector([0 for i in allpos]),
                          'name': robjects.StrVector([reg_name for i in allpos])}
            for suffix,name in strands.iteritems():
                connection = sqlite3.connect(db+suffix)
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
#            for n in range(prod_shift,len(allpos)-prod_shift):
#                data_block['prod'][n] = sqrt(data_block['plus'][n-prod_shift]*
#                                             data_block['minus'][n+prod_shift]) 
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

