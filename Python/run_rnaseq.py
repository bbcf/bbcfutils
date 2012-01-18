#!/usr/bin/env python
"""
A High-throughput RNA-seq analysis workflow.

python run_rnaseq.py -v lsf -c config_files/gapdh.txt -d rnaseq -p transcripts
python run_rnaseq.py -v lsf -c config_files/rnaseq.txt -d rnaseq -p genes -m /scratch/cluster/monthly/jdelafon/mapseq
"""
import os, sys, json, re
import optparse
from bbcflib import rnaseq, frontend, common, mapseq, genrep, email, gdv
from bein.util import use_pickle, add_pickle
from bein import execution, MiniLIMS

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main():
    map_args = None # {'bwt_args':["-n",str(3),"-p",str(4),"-d",str(50),"--chunkmbs",str(1024),"-m",str(5)]}

    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key of the new RNA-seq job", {'default': None}),
            ("-d", "--rnaseq-minilims", "MiniLIMS where RNAseq executions and files will be stored.",
                                     {'default': None, 'dest':"minilims"}),
            ("-m", "--mapseq-minilims", "MiniLIMS where a previous Mapseq execution and files has been stored. \
                                     Set it to None to align de novo from read files.",
                                     {'default': "/data/htsstation/mapseq/mapseq_minilims", 'dest':"ms_limspath"}),
            ("-w", "--working-directory", "Create execution working directories in wdir",
                                     {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}),
            ("-p", "--pileup_level", "Target features, inside of quotes, separated by commas.\
                                     E.g. 'genes,exons,transcripts'",{'default': "genes,exons,transcripts"}),
            ("-u", "--unmapped", "If True, add junction reads to the pileups.",
                                     {'action':'store_true', 'default': False}))
    try:
        usage = "run_rnaseq.py [OPTIONS]"
        desc = """A High-throughput RNA-seq analysis workflow. It returns a file containing
                  a column of transcript counts for each given BAM file, normalized using DESeq's
                  size factors. """
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (opt, args) = parser.parse_args()

        if os.path.exists(opt.wdir): os.chdir(opt.wdir)
        else: parser.error("Working directory '%s' does not exist." % opt.wdir)
        if not opt.minilims: parser.error("Must specify a MiniLIMS to attach to")

        # Rna-seq job configuration
        M = MiniLIMS(opt.minilims)
        if opt.key:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_rnaseq']['url'] )
            job = htss.job(opt.key) # new *RNA-seq* job instance
            [M.delete_execution(x) for x in M.search_executions(with_description=opt.key,fails=True)]
            description = "Job run with mapseq key %s" % opt.key
        elif os.path.exists(opt.config):
            (job,gl) = frontend.parseConfig(opt.config)
            description = "Job run with config file %s" % opt.config
        else: raise ValueError("Need either a job key (-k) or a configuration file (-c).")
        pileup_level = opt.pileup_level.split(',')

        job.options['unmapped'] = job.options.get('unmapped',True)
        job.options['ucsc_bigwig'] = job.options.get('ucsc_bigwig',True)
        job.options['create_gdv_project'] = job.options.get('create_gdv_project',False)
        if isinstance(job.options['create_gdv_project'],str):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        job.options['discard_pcr_duplicates'] = False
        g_rep = genrep.GenRep( gl.get('genrep_url'), gl.get('bwt_root') )
            #intype is for mapping on the genome (intype=0), exons (intype=1) or transcriptome (intype=2)
        assembly = genrep.Assembly(assembly=job.assembly_id, genrep=g_rep, intype=1)

        # Retrieve mapseq output
        mapseq_url = None
        if 'hts_mapseq' in gl: mapseq_url = gl['hts_mapseq']['url']

        # Program body #
        with execution(M, description=description, remote_working_directory=opt.wdir ) as ex:
            print "Current working directory:", ex.working_directory
            if opt.ms_limspath == "None":
                print "Alignment..."
                job = mapseq.get_fastq_files( job, ex.working_directory)
                fastq_root = os.path.abspath(ex.working_directory)
                bam_files = mapseq.map_groups(ex, job, fastq_root, assembly_or_dict=assembly, map_args=map_args)
                print "Reads aligned."
            else:
                print "Loading BAM files..."
                (bam_files, job) = mapseq.get_bam_wig_files(ex, job, minilims=opt.ms_limspath, hts_url=mapseq_url,
                         script_path=gl.get('script_path') or '', via=opt.via, fetch_unmapped=opt.unmapped)
                assert bam_files, "Bam files not found."
                print "Loaded."
            rnaseq.rnaseq_workflow(ex, job, assembly, bam_files, pileup_level=pileup_level, via=opt.via, unmapped=opt.unmapped)
            gdv_project = {}
            if job.options['create_gdv_project']:
                gdv_project = gdv.new_project( gl['gdv']['email'], gl['gdv']['key'],
                                               job.description, assembly.id, gl['gdv']['url'] )
                add_pickle( ex, gdv_project, description=common.set_file_descr("gdv_json",step='gdv',type='py',view='admin') )
        # GDV
        allfiles = common.get_files(ex.id, M)
        if job.options['create_gdv_project'] and re.search(r'success',gdv_project.get('message','')):
            gdv_project_url = gl['gdv']['url']+"public/project?k="+str(gdv_project['project']['key']) \
                              +"&id="+str(gdv_project['project']['id'])
            allfiles['url'] = {gdv_project_url: 'GDV view'}
            download_url = gl['hts_rnaseq']['download']
            urls  = [download_url+str(k) for k in allfiles['sql'].keys()]
            names = [re.sub('\.sql.*','',str(f)) for f in allfiles['sql'].values()]
            for nurl,url in enumerate(urls):
                try:
                    gdv.new_track( gl['gdv']['email'], gl['gdv']['key'],
                                   project_id=gdv_project['project']['id'],
                                   url=url, file_names=names[nurl],
                                   serv_url=gl['gdv']['url'] )
                except:
                    pass
        print json.dumps(allfiles)

        # E-mail
        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email),
                                   subject="RNA-seq job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''Your RNA-seq job is finished.
                \n The description was: '''+str(job.description)+''' and its unique key is '''+opt.key+'''.
                \n You can retrieve the results at this url: '''+gl['hts_rnaseq']['url']+"jobs/"+opt.key+"/get_results" )
            r.send()

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2


if __name__ == '__main__':
    sys.exit(main())

