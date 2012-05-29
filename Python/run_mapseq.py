#!/usr/bin/env python

"""
A High-throughput sequencing data mapping workflow.

"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle
from bbcflib import daflims, genrep, frontend, email, gdv
from bbcflib.common import get_files, set_file_descr, track_header
from bbcflib.mapseq import *
import sys, optparse, os, re, json

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main():
    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key specifying the job", {'default': None}),
            ("-d", "--minilims", "MiniLIMS where mapseq executions and files will be stored.", {'default': None}),
            ("-w", "--working-directory", "Create execution working directories in wdir", {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}))
    try:
        usage = "run_mapseq.py [-h] [-u via] [-w wdir] [-k job_key] [-c config_file] -d minilims"
        desc = """A High-throughput sequencing data mapping workflow."""
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
            htss = frontend.Frontend( url=gl['hts_mapseq']['url'] )
            job = htss.job( opt.key )
            [M.delete_execution(x) for x in M.search_executions(with_description=opt.key,fails=True)]
            if opt.config and os.path.exists(opt.config):
                (job,gl) = frontend.parseConfig( opt.config, job, gl )
        elif os.path.exists(opt.config):
            (job,gl) = frontend.parseConfig( opt.config )
            opt.key = job.description
        else:
            raise Usage("Need either a job key (-k) or a configuration file (-c).")
        g_rep = genrep.GenRep( url=gl.get("genrep_url"), root=gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep,
                                    intype=job.options.get('input_type_id',0) )
        if 'lims' in gl:
            dafl = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'], password=pwd ))
                        for loc,pwd in gl['lims']['passwd'].iteritems())
        else:
            dafl = None
        job.options['compute_densities'] = job.options.get('compute_densities',True)
        if isinstance(job.options['compute_densities'],basestring):
            job.options['compute_densities'] = job.options['compute_densities'].lower() in ['1','true','t']
        job.options['ucsc_bigwig'] = job.options.get('ucsc_bigwig',True)
        if isinstance(job.options['ucsc_bigwig'],basestring):
            job.options['ucsc_bigwig'] = job.options['ucsc_bigwig'].lower() in ['1','true','t']
        job.options['ucsc_bigwig'] &= (assembly.intype == 0)
        job.options['create_gdv_project'] = job.options.get('create_gdv_project',False)
        if isinstance(job.options['create_gdv_project'],basestring):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        gdv_project = {'project':{'id': job.options.get('gdv_project_id',0)}}
        if not('gdv_key' in job.options): job.options['gdv_key'] = ""
        job.options['create_gdv_project'] = job.options['create_gdv_project'] and job.options['compute_densities']
        map_args = job.options.get('map_args',{})
        map_args['via'] = opt.via
        logfile = open(opt.key+".log",'w')
        debugfile = open(opt.key+".debug",'w')
        debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()
        with execution( M, description=opt.key, remote_working_directory=opt.wdir ) as ex:
            logfile.write("Enter execution, fetch fastq files.\n");logfile.flush()
            job = get_fastq_files( ex, job, dafl )
            logfile.write("Generate QC report.\n");logfile.flush()
            run_fastqc( ex, job, via=opt.via )
            logfile.write("Map reads.\n");logfile.flush()
            mapped_files = map_groups( ex, job, assembly, map_args )
            logfile.write("Make stats:\n");logfile.flush()
            debugfile.write("GroupId_GroupName:\t")
            for k,v in job.groups.iteritems():
                debugfile.write(str(k)+"_"+str(v['name'])+"\t");debugfile.flush()
                pdf = add_pdf_stats( ex, mapped_files,
                                     {k:v['name']},
                                     gl.get('script_path',''),
                                     description=set_file_descr(v['name']+"_mapping_report.pdf",groupId=k,step='stats',type='pdf') )
            if job.options['compute_densities']:
                logfile.write("\ncomputing densities.\n");logfile.flush()
                if int(job.options.get('read_extension',-1))<=0:
                    job.options['read_extension'] = mapped_files.values()[0].values()[0]['stats']['read_length']
                density_files = densities_groups( ex, job, mapped_files, assembly.chromosomes, via=opt.via )
                logfile.write("Finished computing densities.\n");logfile.flush()
                if job.options['create_gdv_project']:
                    gdv_project = gdv.get_project(mail=gl['gdv']['email'], key=gl['gdv']['key'], project_key=job.options['gdv_key'])
                    if 'error' in gdv_project:
                        logfile.write("Creating GDV project.\n");logfile.flush()
                        gdv_project = gdv.new_project( gl['gdv']['email'], gl['gdv']['key'],
                                                       job.description, assembly.id, gl['gdv']['url'] )
                    debugfile.write("\nGDV project: "+json.dumps(gdv_project)+"\n");debugfile.flush()
                    add_pickle( ex, gdv_project, description=set_file_descr("gdv_json",step='gdv',type='py',view='admin') )
        allfiles = get_files( ex.id, M )
        if job.options['ucsc_bigwig']:
            logfile.write("UCSC track file: "+opt.key+".bed\n");logfile.flush()
            ucscfiles = get_files( ex.id, M, select_param={'ucsc':'1'} )
            with open(opt.key+".bed",'w') as ucscbed:
                for ftype,fset in ucscfiles.iteritems():
                    for ffile,descr in fset.iteritems():
                        if re.search(r' \(.*\)',descr): continue
                        ucscbed.write(track_header(descr,ftype,gl['hts_mapseq']['download'],ffile))
        if gdv_project.get('project',{}).get('id',0)>0:
            gdv_project_url = gl['gdv']['url']+"public/project?k="+str(gdv_project['project']['download_key'])+"&id="+str(gdv_project['project']['id'])
            allfiles['url'] = {gdv_project_url: 'GDV view'}
            download_url = gl['hts_mapseq']['download']
            urls  = [download_url+str(k) for k in allfiles['sql'].keys()]
            names = [re.sub('\.sql.*','',str(f)) for f in allfiles['sql'].values()]
            logfile.write("Uploading GDV tracks:\n"+" ".join(urls)+"\n"+" ".join(names)+"\n");logfile.flush()
            tr = gdv.multiple_tracks(mail=gl['gdv']['email'], key=gl['gdv']['key'], serv_url=gl['gdv']['url'], 
                                     project_id=gdv_project['project']['id'], 
                                     extensions=['sql']*len(urls),
                                     urls=urls, tracknames=names, force=True )
            debugfile.write("GDV Tracks Status\n"+"\n".join([str(v) for v in tr])+"\n");debugfile.flush()
        logfile.close()
        debugfile.close()
        print json.dumps(allfiles)
        with open(opt.key+".done",'w') as done:
            json.dump(allfiles,done)
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email),
                                   subject="Mapseq job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''
Your mapseq job has finished.

The description was:
'''+str(job.description)+'''
and its unique key is '''+opt.key+'''.

You can now retrieve the results at this url:
'''+gl['hts_mapseq']['url']+"jobs/"+opt.key+"/get_results")
            r.send()
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
