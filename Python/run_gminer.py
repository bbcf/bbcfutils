#!/usr/bin/env python

import sys
import pickle
import gMiner
import optparse


class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
           )
    try:
        usage = "run_gminer.py pickle_file [OPTIONS]"
        desc = """ Runs gMiner [...]
               Arguments:
               `pickle_file`: a pickle file with a dictionary describing the gMiner job,
               saves the output as an entry 'job_output' in the same dictionary, same file. """
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (opt, args) = parser.parse_args()

        if len(args) != 1:
            parser.error("run_gminer.py takes exactly 1 argument.")
        pickle_file = args[0]
        with open(pickle_file,'r') as f:
            job = pickle.load(f)
        job['job_output'] = gMiner.run(**job)
        with open(pickle_file,'w') as f:
            pickle.dump(job,f)
        sys.exit(0)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())
