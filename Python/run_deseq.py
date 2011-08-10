#!/bin/env python
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

usage = """run_deseq.py cond1_file cond2_file cond1_label cond2_label method [output] [gene_names]

*cond1_file*, *cond2_file*:     arrays of counts
*cond1_label*, *cond2_label* :  string describing the respective condition
*method*:                       'normal' or 'blind'
*output*:                       name of the output file (optional)
*gene_names*:                   dictionary {gene ID: gene name} (optional)
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    """
    Usage: run_deseq.py cond1_file cond2_file cond1_label cond2_label method [output] [gene_names]
    """
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) < 5:
            raise Usage("run_deseq.py takes at least %d arguments." % 5)
        cond1_file = argv[0]
        with open(cond1_file,'rb') as f:
            cond1 = pickle.load(f)
        cond2_file = argv[1]
        with open(cond2_file,'rb') as f:
            cond2 = pickle.load(f)
        cond1_label = str(argv[2])
        cond2_label = str(argv[3])
        method = str(argv[4])
        output = str(argv[5])
        gene_names_file = argv[6]
        with open(gene_names_file,'rb') as f:
            gene_names = pickle.load(f)

        optargs = {}
        if not output=="None": optargs['output'] = output
        if not gene_names=="None": optargs['gene_names'] = gene_names
        
        result_filename = inference(cond1_label, cond1, cond2_label, cond2, method, **optargs)
        
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())


