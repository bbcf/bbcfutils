#!/bin/env python
"""
A High-throughput sequencing data mapping workflow.

"""
from bbcflib import daflims, genrep, frontend, email, gdv, common
from bbcflib.mapseq import *
import sys, getopt, os

usage = """run_mapseq.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] -d minilims

-h           Print this message and exit
-u via       Run executions using method 'via' (can be "local" or "lsf")
-w wdir      Create execution working directories in wdir
-d minilims  MiniLIMS where mapseq executions and files will be stored.
-k job_key   Alphanumeric key specifying the job
-c file      Config file
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    via = "lsf"
    limspath = None
    hts_key = None
    working_dir = None
    config = None
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hu:k:d:w:",
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
        M = MiniLIMS( limspath )
        if len(hts_key)>1:
            gl = use_pickle(M, "global variables")
            htss = frontend.Frontend( url=gl['hts_mapseq']['url'] )
            job = htss.job( hts_key )
        ###[M.delete_execution(x) for x in M.search_executions(with_description=hts_key)]
        elif os.path.exists(config_file):
            (job,gl) = frontend.parseConfig( config_file )
        else:
            raise ValueError("Need either a job key (-k) or a configuration file (-c).")
        g_rep = genrep.GenRep( gl["genrep_url"], gl["bwt_root"] )
        assembly = g_rep.assembly( job.assembly_id )
        dafl = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'], password=pwd ))
                    for loc,pwd in gl['lims']['passwd'].iteritems())
        job.options['ucsc_bigwig'] = job.options.get('ucsc_bigwig') or True
        job.options['gdv_project'] = job.options.get('gdv_project') or False
        with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
            job = get_fastq_files( job, ex.working_directory, dafl )
            mapped_files = map_groups( ex, job, ex.working_directory, assembly, {'via': via} )
            pdf = add_pdf_stats( ex, mapped_files,
                                 dict((k,v['name']) for k,v in job.groups.iteritems()),
                                 gl['script_path'] )
            if job.options['compute_densities']:
                if not(job.options.get('read_extension')>0):
                    job.options['read_extension'] = mapped_files.values()[0].values()[0]['stats']['read_length']
                density_files = densities_groups( ex, job, mapped_files, assembly.chromosomes, via=via )
                if job.options['gdv_project']:
                    gdv_project = gdv.create_gdv_project( gl['gdv']['key'], gl['gdv']['email'],
                                                          job.description, hts_key, 
                                                          assembly.nr_assembly_id,
                                                          gdv_url=gl['gdv']['url'], public=True )
                    add_pickle( ex, gdv_project, description='py:gdv_json' )
        allfiles = common.get_files( ex.id, M )
        if 'py:gdv_json' in allfiles:
            allfiles['url'] = {gdv_project['public_url']: 'GDV view'}
            download_url = gl['hts_mapseq']['download']
            [gdv.add_gdv_sqlite( gl['gdv']['key'], gl['gdv']['email'],
                                 gdv_project['project_id'],
                                 url=download_url+str(k), 
                                 name = re.sub('\.sql','',str(f)),
                                 gdv_url=gl['gdv']['url'], datatype="quantitative" ) 
             for k,f in allfiles['sql'].iteritems()]
        print json.dumps(allfiles)
        r = email.EmailReport( sender=gl['email']['sender'],
                               to=str(job.email),
                               subject="Mapseq job "+str(job.description),
                               smtp_server=gl['email']['smtp'] )
        r.appendBody('''
Your mapseq job is finished.

The description was: 
'''+str(job.description)+'''
and its unique key is '''+hts_key+'''.

You can now retrieve the results at this url:
'''+gl['hts_mapseq']['url']+"jobs/"+hts_key+"/get_results")
        r.send()
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    
if __name__ == '__main__':
    sys.exit(main())
