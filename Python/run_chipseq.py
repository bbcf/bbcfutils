#!/usr/bin/env python

"""
A High-throughput ChIP-seq peak analysis workflow.

"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle
from bbcflib import genrep, frontend, email, gdv, mapseq
from bbcflib.common import get_files
from bbcflib.chipseq import *
import sys, optparse, os, json, re

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key specifying the job", {'default': None}),
            ("-d", "--minilims", "MiniLIMS where chipseq executions and files will be stored.", {'default': None}),
            ("-m", "--mapseq_minilims", "MiniLIMS where a previous Mapseq execution and files has been stored.",
             {'default': "/srv/mapseq/public/data/mapseq_minilims"}),
            ("-w", "--working-directory", "Create execution working directories in wdir", {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}))
    try:
        usage = "run_chipseq.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] [-m mapseq_minilims] -d minilims"
        desc = """A ChIP-seq peak analysis workflow."""
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (opt, args) = parser.parse_args()

        if os.path.exists(opt.wdir): os.chdir(opt.wdir)
        else: parser.error("Working directory '%s' does not exist." % opt.wdir)
        if not(opt.minilims and os.path.exists(opt.minilims)
               and (opt.key != None or (opt.config and os.path.exists(opt.config)))):
            raise Usage("Need a minilims and a job key or a configuration file")
        M = MiniLIMS( opt.minilims )
        if opt.key:
            gl = use_pickle(M, "global variables")
            htss = frontend.Frontend( url=gl['hts_chipseq']['url'] )
            job = htss.job( opt.key )
            [M.delete_execution(x) for x in M.search_executions(with_description=opt.key,fails=True)]
            if opt.config and os.path.exists(opt.config):
                (job,gl) = frontend.parseConfig( opt.config, job, gl )
        elif os.path.exists(opt.config):
            (job,gl) = frontend.parseConfig( opt.config )
            opt.key = job.description
        else:
            raise Usage("Need either a job key (-k) or a configuration file (-c).")
        mapseq_url = None
        if 'hts_mapseq' in gl:
            mapseq_url = gl['hts_mapseq']['url']

        job.options.setdefault('ucsc_bigwig',True)
        if isinstance(job.options['ucsc_bigwig'],basestring):
            job.options['ucsc_bigwig'] = job.options['ucsc_bigwig'].lower() in ['1','true','t']
        job.options.setdefault('create_gdv_project',False)
        if isinstance(job.options['create_gdv_project'],basestring):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        gdv_project = {'project':{'id': job.options.get('gdv_project_id',0)}}
        job.options.setdefault('gdv_key',"")

        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )
        logfile = open(opt.key+".log",'w')
        debugfile = open(opt.key+".debug",'w')
        debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()
        with execution( M, description=opt.key, remote_working_directory=opt.wdir ) as ex:
            logfile.write("Enter execution, fetch bam and wig files.\n");logfile.flush()
            job = mapseq.get_bam_wig_files( ex, job, minilims=opt.mapseq_minilims, hts_url=mapseq_url,
                                                            script_path=gl.get('script_path',''), via=opt.via )
            logfile.write("Starting workflow.\n");logfile.flush()
            chipseq_files = workflow_groups( ex, job, assembly,
                                             gl.get('script_path',''), logfile=logfile, via=opt.via )
            if job.options.get('create_gdv_project'):
                gdv_project = gdv.get_project(mail=gl['gdv']['email'], key=gl['gdv']['key'], project_key=job.options['gdv_key'])
                if 'error' in gdv_project:
                    logfile.write("\nCreating GDV project.\n");logfile.flush()
                    gdv_project = gdv.new_project( gl['gdv']['email'], gl['gdv']['key'],
                                                   job.description, assembly.id,
                                                   gl['gdv']['url'] )
                debugfile.write("GDV project: "+str(gdv_project['project']['id'])+"\n");debugfile.flush()
                add_pickle( ex, gdv_project, description=set_file_descr("gdv_json",step='gdv',type='py',view='admin') )
        allfiles = get_files( ex.id, M )
        if gdv_project.get('project',{}).get('id',0)>0 and 'sql' in allfiles:
            gdv_project_url = gl['gdv']['url']+"public/project?k="+str(gdv_project['project']['download_key'])+"&id="+str(gdv_project['project']['id'])
            allfiles['url'] = {gdv_project_url: 'GDV view'}
            download_url = gl['hts_chipseq']['download']
            urls  = [download_url+str(k) for k in allfiles['sql'].keys()]
            names = [re.sub('\.sql.*','',str(f)) for f in allfiles['sql'].values()]
            logfile.write("Uploading GDV tracks:\n"+" ".join(urls)+"\n"+" ".join(names)+"\n");logfile.flush()
            try:
                tr = gdv.multiple_tracks(mail=gl['gdv']['email'], key=gl['gdv']['key'], serv_url=gl['gdv']['url'],
                                         project_id=gdv_project['project']['id'],
                                         extensions=['sql']*len(urls),
                                         urls=urls, tracknames=names, force=True )
                debugfile.write("GDV Tracks Status\n"+"\n".join([str(v) for v in tr])+"\n");debugfile.flush()
            except Exception, err:
                debugfile.write("GDV Tracks Failed: %s\n" %err);debugfile.flush()
                pass
        logfile.close()
        debugfile.close()
        print json.dumps(allfiles)
        with open(opt.key+".done",'w') as done:
            json.dump(allfiles,done)
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email).split(','),
                                   subject="Chipseq job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''
Your chip-seq job has finished.

The description was:
'''+str(job.description)+'''
and its unique key is '''+opt.key+'''.

You can retrieve the results at this url:
'''+gl['hts_chipseq']['url']+"jobs/"+opt.key+"/get_results")
            r.send()
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
