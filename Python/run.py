#!/usr/bin/env python

"""
Job launcher.
"""
import os, sys, json, re, sqlite3
import optparse
from bbcflib import frontend, common, email, gdv
from bein.util import use_pickle, add_pickle
from bein import execution, MiniLIMS

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def run(workflow_object):
    WO = workflow_object
    try:
        parser = optparse.OptionParser(usage=WO.usage, description=WO.desc)
        for opt in WO.opts:
            if len(opt)==4:
                parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
            elif len(opt)==3:
                parser.add_option(opt[0],help=opt[1],**opt[2])
        (opt, args) = parser.parse_args()

        if os.path.exists(opt.wdir): os.chdir(opt.wdir)
        else: parser.error("Working directory '%s' does not exist." % opt.wdir)
        if not(opt.minilims and os.path.exists(opt.minilims)):
            parser.error("Need to specify an existing minilims (-d), got %s." % opt.minilims)

        # Job configuration #
        try: M = MiniLIMS(opt.minilims)
        except sqlite3.OperationalError:
            raise ValueError("MiniLIMS not found: %s ." % os.path.abspath(opt.minilims))
        if opt.key:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl[WO.hts]['url'] )
            job = htss.job(opt.key) # new Job instance
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

        logfile = open((opt.key or opt.config)+".log",'w')
        debugfile = open((opt.key or opt.config)+".debug",'w')
        debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()

        # Program body #
        with execution(M, description=description, remote_working_directory=opt.wdir ) as ex:
            logfile.write("Enter execution. Current working directory: %s \n" % ex.working_directory);logfile.flush()

            ### MAIN #############################
            WO.workflow(ex,job,opt,gl,logfile,debugfile)
            ######################################

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
            download_url = gl[WO.hts]['download']
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
                                   subject=WO.name + str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''Your %s is finished.
                \n The description is: %s and its unique key is: %s .
                \n You can retrieve the results at this url: %s''' \
                % (WO.name,job.description,opt.key,gl[WO.hts]['url']+"jobs/"+opt.key+"/get_results"))
            r.send()

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, WO.usage
        return 2

#------------------------------------------------------#
# This code was written by                             #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#
