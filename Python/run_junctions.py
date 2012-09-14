#!/usr/bin/env python

"""
A de novo splice junctions detection workflow.

python run_junctions.py -v lsf -c config_files/junc.config -d junctions
"""
import os, sys, json, re
import optparse
from bbcflib import junctions, frontend, common, mapseq, email, gdv
from bein.util import use_pickle, add_pickle
from bein import execution, MiniLIMS


class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main():
    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key of the new RNA-seq job", {'default': None}),
            ("-d", "--junctions_minilims", "MiniLIMS where RNAseq executions and files will be stored.",
                                     {'default': "/srv/junctions/public/data/junctions_minilims"}),
            ("-m", "--mapseq_minilims", "MiniLIMS where a previous Mapseq execution and files has been stored.",
                                     {'default': "/srv/mapseq/public/data/mapseq_minilims"}),
            ("-w", "--working-directory", "Create execution working directories in wdir",
                                     {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}),
           )
    try:
        usage = "run_junctions.py [-h -v via {-k key | -c config_file} -w working_directory -d minilims -m mapseq_minilims]"
        desc = """Uses SOAPsplice to find transcriptional splice junctions."""
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            if len(opt)==4:
                parser.add_option(opt[0],opt[1],help=opt[2],**opt[3]) # -x / --xxx option
            elif len(opt)==3:
                parser.add_option(opt[0],help=opt[1],**opt[2])        # only --xxx option
        (opt, args) = parser.parse_args()

        if os.path.exists(opt.wdir): os.chdir(opt.wdir)
        else: parser.error("Working directory '%s' does not exist." % opt.wdir)
        if not(opt.junctions_minilims and os.path.exists(opt.junctions_minilims)
                and (opt.key != None or (opt.config and os.path.exists(opt.config)))):
            parser.error("Need a minilims and either a job key (-k) or a configuration file (-c).\n")

        # Junctions job configuration #
        M = MiniLIMS(opt.junctions_minilims)
        if opt.key:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_junctions']['url'] )
            job = htss.job(opt.key) # new *RNA-seq* job instance
            [M.delete_execution(x) for x in M.search_executions(with_description=opt.key,fails=True)]
        elif opt.config and os.path.exists(opt.config):
            (job,gl) = frontend.parseConfig(opt.config)
        else: raise Usage("Need either a job key (-k) or a configuration file (-c).")
        description = opt.key or opt.config

        job.options.setdefault('ucsc_bigwig',True)
        job.options.setdefault('create_gdv_project',False)
        job.options.setdefault('gdv_key',"")
        if isinstance(job.options['create_gdv_project'],str):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        gdv_project = {'project':{'id': job.options.get('gdv_project_id',0)}}
        job.options['discard_pcr_duplicates'] = False

        logfile = open((opt.key or opt.config)+".log",'w')
        debugfile = open((opt.key or opt.config)+".debug",'w')
        debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()

        # Retrieve mapseq output
        mapseq_url = gl.get('hts_mapseq',{}).get('url')

        # Program body #
        with execution(M, description=description, remote_working_directory=opt.wdir ) as ex:
            debugfile.write("Enter execution. Current working directory: %s \n" % ex.working_directory);debugfile.flush()
            logfile.write("Fetch bam and wig files\n");logfile.flush()
            (bam_files, job) = mapseq.get_bam_wig_files(ex, job, minilims=opt.mapseq_minilims, hts_url=mapseq_url,
                                     script_path=gl.get('script_path',''), via=opt.via, fetch_unmapped=True)
            assert bam_files, "Bam files not found."
            logfile.write("Starting workflow.\n");logfile.flush()
            junctions.junctions_workflow(ex, job, bam_files, via=opt.via)

            # Create GDV project #
            if job.options['create_gdv_project']:
                gdv_project = gdv.get_project(mail=gl['gdv']['email'], key=gl['gdv']['key'],
                                              project_key=job.options['gdv_key'])
                if 'error' in gdv_project:
                    logfile.write("Creating GDV project.\n");logfile.flush()
                    gdv_project = gdv.new_project( gl['gdv']['email'], gl['gdv']['key'],
                                                   job.description, job.assembly_id, gl['gdv']['url'] )
                debugfile.write("GDV project: "+json.dumps(gdv_project)+"\n");debugfile.flush()
                add_pickle(ex, gdv_project, description=common.set_file_descr("gdv_json",step='gdv',type='py',view='admin'))

        # Upload tracks to GDV #
        allfiles = common.get_files(ex.id, M)
        if gdv_project.get('project',{}).get('id',0)>0:
            gdv_project_url = gl['gdv']['url']+"public/project?k="+str(gdv_project['project']['download_key']) \
                              +"&id="+str(gdv_project['project']['id'])
            allfiles['url'] = {gdv_project_url: 'GDV view'}
            download_url = gl['hts_junctions']['download']
            urls  = [download_url+str(k) for k in allfiles['sql'].keys()]
            names = [re.sub('\.sql.*','',str(f)) for f in allfiles['sql'].values()]
            logfile.write("Uploading GDV tracks:\n"+" ".join(urls)+"\n"+" ".join(names)+"\n");logfile.flush()
            try:
                tr = gdv.multiple_tracks(mail=gl['gdv']['email'], key=gl['gdv']['key'], serv_url=gl['gdv']['url'],
                                         project_id=gdv_project['project']['id'],
                                         urls=urls, tracknames=names, force=True )
                debugfile.write("GDV Tracks Status\n"+"\n".join([str(v) for v in tr])+"\n");debugfile.flush()
            except Exception, err:
                debugfile.write("GDV Tracks Failed: %s\n" %err);debugfile.flush()
                pass
        logfile.close()
        debugfile.close()
        print json.dumps(allfiles)
        with open((opt.key or opt.config)+".done",'w') as done:
            json.dump(allfiles,done)

        # E-mail #
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email).split(','),
                                   subject="Junctions job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''Your Junctions job is finished.
                \n The description was: '''+str(job.description)+''' and its unique key is '''+opt.key+'''.
                \n You can retrieve the results at this url: '''+gl['hts_junctions']['url']+"jobs/"+opt.key+"/get_results" )
            r.send()

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2


if __name__ == '__main__':
    sys.exit(main())


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#
