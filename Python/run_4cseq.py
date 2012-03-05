#!/usr/bin/env python

"""
A High-throughput 4C-seq analysis workflow.

"""
from bbcflib import genrep, frontend, email, gdv, common, mapseq
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle
import sys, getopt, os, re, json

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
            [M.delete_execution(x) for x in M.search_executions(with_description=hts_key,fails=True)]
        elif os.path.exists(config_file):
            (job,gl) = frontend.parseConfig( config_file )
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
        gdv_project = {'project':{'id': job.options.get('gdv_project_id',0)}}
        if not('gdv_key' in job.options): job.options['gdv_key'] = ""
        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )
        primers_file=os.path.join(working_dir,'primers.fa')
        primers_dict=c4seq.loadPrimers(primers_file)
	logfile = open(hts_key+".log",'w')
	debugfile = open(hts_key+".debug",'w')
	debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()
        with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
	    logfile.write("Enter execution, fetch bam and wig files.\n");logfile.flush()
            (mapseq_files, job) = mapseq.get_bam_wig_files( ex, job, ms_limspath, mapseq_url, suffix=['merged'],script_path=gl['script_path'], via=via )
            logfile.write("Starting workflow.\n");logfile.flush()
	    c4seq_files = c4seq.workflow_groups( ex, job, primers_dict, assembly,
                                                 mapseq_files, mapseq_url,
                                                 gl['script_path'],logfile=logfile)
            if job.options.get('create_gdv_project'):
		gdv_project=gdv.get_project(mail=gl['gdv']['email'], key=gl['gdv']['key'], project_key=job.options['gdv_key'])
		if 'error' in gdv_project:
                    logfile.write("Creating GDV project.\n");logfile.flush()
                    gdv_project = gdv.new_project( gl['gdv']['email'], gl['gdv']['key'],
                                                   job.description, assembly.id, 
                                                   gl['gdv']['url'] )
                debugfile.write("GDV project: "+str(gdv_project['project']['id'])+"\n");debugfile.flush()
                add_pickle( ex, gdv_project, description=common.set_file_descr("gdv_json",step='gdv',type='py',view='admin') )

        ucscfiles = common.get_files( ex.id, M, select_param={'ucsc':'1'} )
        with open(hts_key+".bed",'w') as ucscbed:
            for ftype,fset in ucscfiles.iteritems():
                for ffile,descr in fset.iteritems():
                    ucscbed.write(common.track_header(descr,ftype,gl['hts_4cseq']['download'],ffile))

        allfiles = common.get_files( ex.id, M )
        if gdv_project.get('project',{}).get('id',0)>0:
            gdv_project_url = gl['gdv']['url']+"public/project?k="+str(gdv_project['project']['key'])+"&id="+str(gdv_project['project']['id'])
            allfiles['url'] = {gdv_project_url: 'GDV view'}
            download_url = gl['hts_4cseq']['download']
	    urls=[]
	    names=[]
            exts=[]
	    for l,t in allfiles.iteritems():
		    for k,v in allfiles[l].iteritems():
			if re.search(r'gdv:1',v):
				urls.append(download_url+str(k))
				if re.search(r'\.sql',str(v)):
                                    names.append(re.sub('\.sql.*','',str(v)))
                                    exts.append('sql')
				if re.search(r'\.bedGraph',str(v)):
                                    names.append(re.sub('\.bedGraph.*','',str(v)))
                                    exts.append('bedGraph')
	
	    logfile.write("Uploading GDV tracks:\n");logfile.flush()
	    debugfile.write("Uploading GDV tracks:\n"+" ".join(urls)+"\n"+" ".join(names)+"\n");debugfile.flush()
            gdv.multiple_tracks(mail=gl['gdv']['email'], key=gl['gdv']['key'], serv_url=gl['gdv']['url'], 
                                project_id=gdv_project['project']['id'], 
                                urls=urls, tracknames=names, extensions=exts, force=True )
	logfile.close()
	debugfile.close()
        print json.dumps(allfiles)
        with open(hts_key+".done",'w') as done:
            json.dump(allfiles,done)


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

