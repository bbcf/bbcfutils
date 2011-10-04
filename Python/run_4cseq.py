#!/bin/env python
"""
A High-throughput 4C-seq analysis workflow.

"""
from bbcflib import daflims, genrep, frontend, email, gdv, common
from bbcflib.mapseq import *
import sys, getopt, os, json, re

from bbcflib import c4seq

usage = """run_4cseq.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] -d minilims

-h           Print this message and exit
-u via       Run executions using method 'via' (can be "local" or "lsf")
-w wdir      Create execution working directories in wdir
-d minilims  MiniLIMS where 4cseq executions and files will be stored.
-m minilims  MiniLIMS where a previous Mapseq execution and files has been stored.
-k job_key   Alphanumeric key specifying the job
-c file      Config file 
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    via = "lsf"
    limspath = None
    ms_limspath = "/data/htsstation/mapseq/mapseq_minilims"
    hts_key = None
    working_dir = None
    config = None
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
        M = MiniLIMS( limspath )
        if len(hts_key)>1:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_4cseq']['url'] )
            job = htss.job( hts_key )
        elif os.path.exists(config_file):
            (job,gl) = frontend.parseConfig( config_file )
        else:
            raise ValueError("Need either a job key (-k) or a configuration file (-c).")
	mapseq_url = None
        if 'hts_mapseq' in gl:
        	mapseq_url = gl['hts_mapseq']['url']
        job.options['ucsc_bigwig'] = True
	g_rep = genrep.GenRep( gl["genrep_url"], gl.get("bwt_root") )
	assembly = g_rep.assembly( job.assembly_id )
	primers_file='/scratch/cluster/monthly/htsstation/4cseq/'+str(job.id)+'/primers.fa'
	primers_dict=c4seq.loadPrimers(primers_file)
        with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
            (mapseq_files, job) = get_bam_wig_files( ex, job, ms_limspath, mapseq_url, suffix=['merged'],script_path=gl['script_path'], via=via )
	    c4seq_files = c4seq.workflow_groups( ex, job, primers_dict, g_rep,
                                           mapseq_files, mapseq_url,  
                                           gl['script_path'])
        allfiles = common.get_files( ex.id, M )
        job.options['gdv_project'] = job.options.get('gdv_project') or True
        if 'gdv_project' in job.options and 'sql' in allfiles:
            allfiles['url'] = {job.options['gdv_project']['public_url']: 'GDV view'}
            download_url = gl['hts_4cseq']['download']
            if job.options['gdv_project']:
                    gdv_project = gdv.create_gdv_project( gl['gdv']['key'], gl['gdv']['email'],
                                                          job.description,
                                                          assembly.nr_assembly_id,
                                                          gdv_url=gl['gdv']['url'], public=True )
                    add_pickle( ex, gdv_project, description='py:gdv_json' )
                    
            [gdv.add_gdv_track( gl['gdv']['key'], gl['gdv']['email'],
                                 job.options['gdv_project']['project_id'],
                                 url=download_url+str(k), 
                                 name = re.sub('\.sql','',str(f)),
                                 gdv_url=gl['gdv']['url']) 
             for k,f in allfiles['sql'].iteritems()]
        print json.dumps(allfiles)

        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email),
                                   subject="4cseq job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''
Your 4C-seq job is finished.

The description was: 
'''+str(job.description)+'''
and its unique key is '''+hts_key+'''.

You can retrieve the results at this url:
'''+gl['hts_4cseq']['url']+"jobs/"+hts_key+"/get_results")
            r.send()
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    
if __name__ == '__main__':
    sys.exit(main())
