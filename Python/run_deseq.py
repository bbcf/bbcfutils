#!/bin/env python
import os, sys, getopt
import pickle
from bbcflib.rnaseq import comparisons

usage = """run_deseq.py cond1_file cond2_file cond1_label cond2_label assembly_id target method output

*cond1_file*, *cond2_file*      arrays (.pickle) of counts
*cond1_label*, *cond2_label*    string describing the respective condition
*assembly_id*                   assembly ID (e.g. 'hg19' or 76)
*target*                        array (.pickle) indicating the target features (e.g. ['genes','exons'])
*method*                        'normal', 'blind' or 'pooled', for DESeq
*output*                        name of the output pickle file
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    """
    Usage: run_deseq.py cond1_file cond2_file cond1_label cond2_label assembly_id target method output
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
        assembly_id = argv[4]
        
        targ = argv[5]
        with open(targ,'rb') as f:
            target = pickle.load(f)
            
        method = str(argv[6])
        output = str(argv[7])

        optargs = {}
        #if not maplot=="None": optargs['maplot'] = maplot
        
        result = comparisons(cond1_label, cond1, cond2_label, cond2, assembly_id, target, method, **optargs)
        with open(output,"wb") as f:
            pickle.dump(result,f)
        
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())


