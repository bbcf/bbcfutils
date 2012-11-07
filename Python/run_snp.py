#!/usr/bin/env python

"""
SNP detection workflow.
"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle
from bbcflib import genrep, frontend, gdv, mapseq, common, snp
from bbcflib import email
import sys, os, json, optparse, re


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv = None):
    opts = (("-v", "--via", "Run executions locally or using bsub (can be either 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key of the new SNP job", {'default': None}),
            ("-m", "--mapseq_limspath", "MiniLIMS where a previous Mapseq execution and files has been stored.",
                                     {'default': "/srv/mapseq/public/data/mapseq_minilims"}),
            ("-w", "--working-directory", "Create execution working directories in wdir",
                                     {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}),
            ("-d", "--snp_limspath", "MiniLIMS where SNP executions and files will be stored.", \
                                     {'default': "/srv/snp/public/data/snp_minilims"}),
            ("-f", "--fasta_path", "Path to a directory containing a fasta file for each chromosome",
                                     {'default':None}),)
    try:
        usage = "run_snp.py [OPTIONS]"
        desc = """Compares sequencing data to a reference assembly to detect SNP."""
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (opt, args) = parser.parse_args()

        if os.path.exists(opt.wdir): os.chdir(opt.wdir)
        else: parser.error("Working directory '%s' does not exist." % opt.wdir)
        if not(os.path.exists(opt.snp_limspath)
               and (opt.key != None or (opt.config and os.path.exists(opt.config)))):
            raise Usage("Need a minilims and a job key or a configuration file")
        M = MiniLIMS( opt.snp_limspath )
        if opt.key:
            hts_key = opt.key
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_snp']['url'] )
            job = htss.job( hts_key )
            [M.delete_execution(x) for x in M.search_executions(with_description=hts_key,fails=True)]
        elif os.path.exists(opt.config):
            (job,gl) = frontend.parseConfig(opt.config)
            hts_key = job.description
        else:
            raise ValueError("Need either a job key (-k) or a configuration file (-c).")
        mapseq_url = gl.get('hts_mapseq',{}).get('url')

        job.options.setdefault('ucsc_bigwig',True)
        job.options.setdefault('create_gdv_project',False)
        job.options.setdefault('gdv_key',"")
        if isinstance(job.options['create_gdv_project'],str):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        gdv_project = {'project':{'id': job.options.get('gdv_project_id',0)}}

        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )

        logfile = open(hts_key+".log",'w')
        debugfile = open(hts_key+".debug",'w')
        debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()

        # Program body
        with execution( M, description=hts_key, remote_working_directory=opt.wdir ) as ex:
            print "Current working directory:", ex.working_directory
            (bam_files, job) = mapseq.get_bam_wig_files(ex, job, minilims=opt.mapseq_limspath,
                                                        hts_url=mapseq_url,
                                                        script_path=gl.get('script_path',''),
                                                        via=opt.via)
            assert bam_files, "Bam files not found."
            logfile.write("cat genome fasta files\n");logfile.flush()
            snp.snp_workflow(ex,job,bam_files,assembly,path_to_ref=opt.fasta_path,via=opt.via)

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
            download_url = gl['hts_snp']['download']
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
        with open(hts_key+".done",'w') as done:
            json.dump(allfiles,done)

        # E-mail #
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email).split(','),
                                   subject="SNP job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''Your SNP job is finished.
                \n The description was: '''+str(job.description)+''' and its unique key is '''+opt.key+'''.
                \n You can retrieve the results at this url: '''+gl['hts_snp']['url']+"jobs/"+opt.key+"/get_results" )
            r.send()

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())



