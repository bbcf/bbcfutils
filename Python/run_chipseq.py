#!/bin/env python
"""
A High-throughput ChIP-seq peak analysis workflow.

"""
from bbcflib import genrep, frontend, email, gdv, common
from bbcflib.mapseq import *
from bbcflib.chipseq import *
import sys, getopt, os, json, re

usage = """run_chipseq.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] -d minilims

-h           Print this message and exit
-u via       Run executions using method 'via' (can be "local" or "lsf")
-w wdir      Create execution working directories in wdir
-d minilims  MiniLIMS where Chipseq executions and files will be stored.
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
    hts_key = ''
    working_dir = None
    config_file = None
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
            elif o in ("-m", "--mapseq_minilims"):
                ms_limspath = a
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
            htss = frontend.Frontend( url=gl['hts_chipseq']['url'] )
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
        g_rep = genrep.GenRep( gl["genrep_url"], gl.get("bwt_root") )
        assembly = g_rep.assembly( job.assembly_id )
        with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
            (mapped_files, job) = get_bam_wig_files( ex, job, minilims=ms_limspath, hts_url=mapseq_url,
                                                     script_path=gl.get('script_path') or '', via=via )
            chipseq_files = workflow_groups( ex, job, mapped_files, assembly.chromosomes, gl.get('script_path') or '' )
            
        allfiles = common.get_files( ex.id, M )
        if 'gdv_project' in job.options and 'sql' in allfiles:
            allfiles['url'] = {job.options['gdv_project']['public_url']: 'GDV view'}
            download_url = gl['hts_chipseq']['download']
            [gdv.add_gdv_track( gl['gdv']['key'], gl['gdv']['email'],
                                job.options['gdv_project']['project_id'],
                                url=download_url+str(k), 
                                name = re.sub('\.sql','',str(f)),
                                gdv_url=gl['gdv']['url'] ) 
             for k,f in allfiles['sql'].iteritems()]
        print json.dumps(allfiles)
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email),
                                   subject="Chipseq job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''
Your chip-seq job is finished.

The description was: 
'''+str(job.description)+'''
and its unique key is '''+hts_key+'''.

You can retrieve the results at this url:
'''+gl['hts_chipseq']['url']+"jobs/"+hts_key+"/get_results")
            r.send()
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    
if __name__ == '__main__':
    sys.exit(main())
