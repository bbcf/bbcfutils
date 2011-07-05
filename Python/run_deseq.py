#!/archive/epfl/bbcf/bin/bin.x86_64/python
import getopt
import os
import sys
import pickle
import numpy
from bbcflib.rnaseq import inference
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc

usage = """run_deseq.py cond1_file cond2_file transcript_names_file cond1_label cond2_label method assembly_id [output] [maplot]

"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    """
    Usage: run_deseq.py cond1_file cond2_file transcript_names_file cond1_label cond2_label method assembly_id [output] [maplot]
    """
    narg = 9
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) < 7:
            raise Usage("run_deseq.py takes at least %d argument." % 7)
        cond1_file = argv[0]
        with open(cond1_file,'rb') as f:
            cond1 = pickle.load(f)
        cond2_file = argv[1]
        with open(cond2_file,'rb') as f:
            cond2 = pickle.load(f)
        transcript_names_file = argv[2]
        with open(transcript_names_file,'rb') as f:
            transcript_names = pickle.load(f)
        cond1_label = str(argv[3])
        cond2_label = str(argv[4])
        method = str(argv[5])
        assembly_id = argv[6]
        output = str(argv[7])
        maplot = str(argv[8])

        optargs = {}
        if not output=="None": optargs['output'] = output
        if not maplot=="None": optargs['maplot'] = maplot
        
        result_filename = inference(cond1_label, cond1, cond2_label, cond2, transcript_names, method, assembly_id, **optargs)
        
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())


