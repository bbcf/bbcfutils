#!/usr/bin/env python

"""
Demultiplexing

"""
from bein import execution, MiniLIMS
from bein.util import use_pickle
from bbcflib import frontend, email, common
from bbcflib.mapseq import *
import sys, getopt, os, json

from bbcflib import demultiplex
#import demultiplex

usage = """run_demultiplex.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] -d minilims

-h           Print this message and exit
-u via       Run executions using method 'via' (can be "local" or "lsf")
-w wdir      Create execution working directories in wdir
-d minilims  MiniLIMS where demultiplexing executions and files will be stored.
-k job_key   Alphanumeric key specifying the job
-c file      Config file
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    via = "lsf"
    limspath = None
    ms_limspath = "/srv/demultiplex/public/data/demultiplex_minilims"
    hts_key = None
    working_dir = os.getcwd()
    config = None
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hu:k:d:w:m:c:",
                                      ["help","via","key","minilims",
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
            elif o in ("-k", "--key"):
                hts_key = a
            elif o in ("-c", "--config"):
                config_file = a
            else:
                raise Usage("Unhandled option: " + o)
    #temporary:  test_demultiplexing.py -i /scratch/cluster/monthly/mleleu/Nicolas/NL1_NoIndex_L004_R1.fastq -p /scratch/cluster/monthly/mleleu/Nicolas/primers4C_GT_NL.fa -x 22 -n 2 -s 77 -l 30 &
    #job={'description':'test_demultiplex_Nico', 'options':{'opt1':'','opt2':''} , 'group': {'grpName':'grpNameNico', 'run':{'runName':'runNameNico','fastaFile':'/scratch/cluster/monthly/mleleu/Nicolas/NL1_NoIndex_L004_R1_part.fastq'}}}
        M = MiniLIMS( limspath )
        if len(hts_key)>1:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_demultiplex']['url'] )
            job = htss.job( hts_key )
            [M.delete_execution(x) for x in M.search_executions(with_description=hts_key,fails=True)]
        elif os.path.exists(config_file):
            (job,gl) = frontend.parseConfig( config_file )
        else:
            raise ValueError("Need either a job key (-k) or a configuration file (-c).")
        job.options['ucsc_bigwig'] = True
        with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
            demultiplex_files = demultiplex.workflow_groups( ex, job, gl)
        allfiles = common.get_files( ex.id, M )

#        gdv_project = gdv.create_gdv_project( gl['gdv']['key'], gl['gdv']['email'],
 #                                               job.description,
  #                                              assembly.nr_assembly_id,
   #                                             gdv_url=gl['gdv']['url'], public=True )
    #    if 'sql' in allfiles:
     #       allfiles['url'] = {gdv_projec['public_url']: 'GDV view'}
     #       download_url = gl['hts_demultiplex']['download']
     #       [gdv.add_gdv_track( gl['gdv']['key'], gl['gdv']['email'],
     #                           gdv_project['project_id'],
     #                           url=download_url+str(k),
     #                           name = re.sub('\.sql','',str(f)),
     #                           gdv_url=gl['gdv']['url'] )
     #        for k,f in allfiles['sql'].iteritems()]
        print json.dumps(allfiles)
        with open(hts_key+".done",'w') as done:
            json.dump(allfiles,done)
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email),
                                   subject="Demultiplexing job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''
Your chip-seq job is finished.

The description was:
'''+str(job.description)+'''
and its unique key is '''+hts_key+'''.

You can retrieve the results at this url:
'''+gl['hts_demultiplex']['url']+"jobs/"+hts_key+"/get_results")
            r.send()
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
