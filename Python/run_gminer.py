#!/archive/epfl/bbcf/bin/bin.x86_64/python
import getopt
import os
import sys
import pickle
import gMiner

usage = """run_gminer.py pickle_file

pickle_file a pickle file with a dictionary describing the gMiner job
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) != 1:
            raise Usage("run_gminer.py takes exactly 1 argument.")
        pickle_file = argv[0]
        with open(pickle_file,'r') as f:
            job = pickle.load(f)
        job['output_files'] = gMiner.run(**job)
        with open(pickle_file,'w') as f:
            pickle.dump(job,f)
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())

