#!/usr/bin/env python

"""
A High-throughput ChIP-seq peak analysis workflow.

"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle
from bbcflib import genrep, frontend, email, gdv, mapseq
from bbcflib.common import get_files
from bbcflib.chipseq import *
import sys, getopt, os, json, re

usage = """run_chipseq.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] [-m mapseq_minilims] -d minilims

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
                                      ["help","via=","key=","minilims=",
                                       "mapseq_minilims=",
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
        if not('create_gdv_project' in job.options):
            job.options['create_gdv_project'] = False
        elif isinstance(job.options['create_gdv_project'],basestring):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )
        logfile = open(hts_key+".log",'w')
        logfile.write(json.dumps(gl)+"\n");logfile.flush()
        with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
            logfile.write("Enter execution, fetch bam and wig files.\n");logfile.flush()
            (mapped_files, job) = mapseq.get_bam_wig_files( ex, job, minilims=ms_limspath, hts_url=mapseq_url,
                                                            script_path=gl.get('script_path') or '', via=via )
            logfile.write("Starting workflow.\n");logfile.flush()
            chipseq_files = workflow_groups( ex, job, mapped_files, assembly,
                                             gl.get('script_path') or '', logfile=logfile, via=via )
            gdv_project = {}
            if job.options.get('create_gdv_project'):
                logfile.write("Creating GDV project.\n");logfile.flush()
                gdv_project = gdv.new_project( gl['gdv']['email'], gl['gdv']['key'], 
                                               job.description, assembly.id,
                                               gl['gdv']['url'] )
                logfile.write("GDV project: "+str(gdv_project['project']['id'])+"\n");logfile.flush()
                add_pickle( ex, gdv_project, description=set_file_descr("gdv_json",step='gdv',type='py',view='admin') )
        allfiles = get_files( ex.id, M )
        if re.search(r'success',gdv_project.get('message','')) and 'sql' in allfiles:
            gdv_project_url = gl['gdv']['url']+"public/project?k="+str(gdv_project['project']['key'])+"&id="+str(gdv_project['project']['id'])
            allfiles['url'] = {gdv_project_url: 'GDV view'}
            download_url = gl['hts_chipseq']['download']
            urls  = [download_url+str(k) for k in allfiles['sql'].keys()]
            names = [re.sub('\.sql.*','',str(f)) for f in allfiles['sql'].values()]
            logfile.write("Uploading GDV tracks:\n"+" ".join(urls)+"\n"+" ".join(names)+"\n");logfile.flush()
            for nurl,url in enumerate(urls):
                try:
                    gdv.new_track( gl['gdv']['email'], gl['gdv']['key'], 
                                   project_id=gdv_project['project']['id'],
                                   url=url, file_names=names[nurl],
                                   serv_url=gl['gdv']['url'] )
                except Exception, e:
                    logfile.write("Error with %s: %s\n" %(names[nurl],e));logfile.flush()
        logfile.close()
        print json.dumps(allfiles)
        with open(hts_key+".done",'w') as done:
            json.dump(allfiles,done)
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email),
                                   subject="Chipseq job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''
Your chip-seq job has finished.

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
