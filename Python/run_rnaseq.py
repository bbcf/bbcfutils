#!/bin/env python
"""
run_rnaseq.py

[FIXME: Description ]
"""
import os
import sys
import getopt
from bbcflib.rnaseq import *
from bbcflib import frontend

usage = """run_rnaseq.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] [-d minilims] [-m minilims]

-h           Print this message and exit
-u via       Run executions using method 'via' (can be "local" or "lsf")
-w wdir      Create execution working directories in wdir
-d minilims  MiniLIMS where RNAseq executions and files will be stored.
-m minilims  MiniLIMS where a previous Mapseq execution and files has been stored.
-k job_key   Alphanumeric key specifying the job
-c file      Config file
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv=None):
    via = "lsf"
    limspath = "rnaseq"
    hts_key = None
    working_dir = os.getcwd()

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hu:k:d:w:m:c:",
                                      ["help","via","key","minilims",
                                       "mapseq_minilims",
                                       "working-directory","config"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-u", "--via"):
                if a=="local":
                    via = "local"
                elif a=="lsf":
                    via = "lsf"
                else:
                    raise Usage("Via (-u) can only be \"local\" or \"lsf\", got %s." % (a,))
            elif o in ("-w", "--working-directory"):
                if os.path.exists(a):
                    os.chdir(a)
                    working_dir = a
                else:
                    raise Usage("Working directory '%s' does not exist." % a)
            elif o in ("-d", "--minilims"):
                limspath = a
            elif o in ("-m", "--mapseq_minilims"):
                ms_limspath = a
            elif o in ("-k", "--key"):
                hts_key = a
            elif o in ("-c", "--config"):
                config_file = a
            else:
                raise Usage("Unhandled option: " + o)

        if len(args) != 0:
            raise Usage("workflow.py takes no arguments without specifiers [-x arg].")
        if limspath == None:
            raise Usage("Must specify a MiniLIMS to attach to")

        M = MiniLIMS(limspath)
        if hts_key:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_rnaseq']['url'] )
            job = htss.job( hts_key )
        elif os.path.exists(config_file):
            (job,gl) = frontend.parseConfig( config_file )
        else:
            raise ValueError("Need either a job key (-k) or a configuration file (-c).")
        job.options['ucsc_bigwig'] = True

        g_rep = genrep.GenRep( gl['genrep_url'], gl['bwt_root'] )
        assembly = g_rep.assembly( assembly_id )

        with execution(M) as ex:
            json = rnaseq_workflow(ex, job, assembly, via=via)

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2


if __name__ == '__main__':
    sys.exit(main())
