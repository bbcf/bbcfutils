#!/usr/bin/env python

"""

SNP detection workflow.

"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle
#from bbcflib import daflims, genrep, frontend, email, gdv
from bbcflib import daflims, genrep, frontend, gdv
from bbcflib.common import get_files, set_file_descr, track_header
from bbcflib.snp import *
import sys, getopt, os, re, json

usage = """run_snp.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] -m [mapseq_minilims] -d minilims

-h Print this message and exit
-u via Run executions using method 'via' (can be "local" or "lsf")
-w wdir Create execution working directories in wdir
-d minilims MiniLIMS where mapseq executions and files will be stored.
-k job_key Alphanumeric key specifying the job
-c file Config file
-m MiniLIMS where a previous Mapseq execution and files has been stored.
"""

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv = None):
    via = "lsf"
    limspath = None
    snp_limspath = "/Users/carat/Desktop/postdoc/bbcf/data/snp_minilims"
    hts_key = ''
    working_dir = None
    config_file = None
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hu:k:d:w:m:c:",
                                      ["help","via=","key=","minilims=",
                                       "snp_minilims=",
                                       "working-directory=","config="])
	except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                return 0
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
            elif o in ("-m", "--snp_minilims"):
                snp_limspath = a
            elif o in ("-k", "--key"):
                hts_key = a
            elif o in ("-c", "--config"):
                config_file = a
            else:
                raise Usage("Unhandled option: " + o)

        if not(limspath and os.path.exists(limspath)
               and (hts_key != None or (config_file and os.path.exists(config_file)))):
            raise Usage("Need a minilims and a job key or a configuration file")
        M = MiniLIMS( limspath )

        if len(hts_key)>1:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_snp']['url'] )
            job = htss.job( hts_key )
            [M.delete_execution(x) for x in M.search_executions(with_description=hts_key,fails=True)]
        elif os.path.exists(config_file):
            (job,gl) = frontend.parseConfig( config_file )
            hts_key = job.description
        else:
            raise ValueError("Need either a job key (-k) or a configuration file (-c).")
        mapseq_url = None

 
        if 'hts_mapseq' in gl:
            mapseq_url = gl['hts_mapseq']['url']
        job.options['ucsc_bigwig'] = True
        if not('create_gdv_project' in job.options):
            job.options['create_gdv_project'] = False
        elif isinstance(job.options['create_gdv_project'],str):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )
        logfile = open(hts_key+".log",'w')
        logfile.write(json.dumps(gl));logfile.flush()

	with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
		logfile.write("test\n")
		#print assembly.id
		(snp_files,job)=snp.pileup(ex,job,minilims=snp_limspath,hts_url=mapseq_url,script_path=gl.get('script_path') or '',via=via)
		logfile.flush()
		logfile.close()
	return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())

