#!/usr/bin/env python

"""

SNP detection workflow.

"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle, pause
from bbcflib import daflims, genrep, frontend, gdv, mapseq, common, snp
from bbcflib import email
from bbcflib.common import set_file_descr
import sys, os, json, optparse


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv = None):
    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key of the new RNA-seq job", {'default': None}),
            ("-m", "--mapseq_limspath", "MiniLIMS where a previous Mapseq execution and files has been stored. \
                                     Set it to None to align de novo from read files.",
                                     {'default': "/data/htsstation/mapseq/mapseq_minilims"}),
            ("-w", "--working-directory", "Create execution working directories in wdir",
                                     {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}),
            ("-d", "--snp_limspath", "MiniLIMS where snp executions and files will be stored.", \
                                     {'default': "/home/scarat/repos/data/snp_minilims"}))
    try:
        usage = "run_snp.py [OPTIONS]"
        desc = """........."""
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
        mapseq_url = None

        if 'hts_mapseq' in gl:
            mapseq_url = gl['hts_mapseq']['url']
        job.options['ucsc_bigwig'] = True
        if not('create_gdv_project' in job.options):
            job.options['create_gdv_project'] = False
        elif isinstance(job.options['create_gdv_project'],str):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )
        logfile = open(hts_key+".log",'w')
        logfile.write(json.dumps(gl)); logfile.flush()

        # Program body
        with execution( M, description=hts_key, remote_working_directory=opt.wdir ) as ex:
            (bam_files, job) = mapseq.get_bam_wig_files(ex, job, minilims=opt.mapseq_limspath, hts_url=mapseq_url, \
                                                        script_path=gl.get('script_path') or '', via=opt.via)
            assert bam_files, "Bam files not found."
            print "cat genome fasta files"
            genomeRef=snp.untar_cat(ex,assembly.fasta_path())
            print "done"

            # for each group
            dictPileupFile={}
            for idGroup,dictRuns in job.groups.iteritems():
                nbRuns=len(dictRuns["runs"].keys())
                # for each run
                listFormatedFile=[]
                for idRun,dictRun in dictRuns["runs"].iteritems():
                    pileupFilename=common.unique_filename_in()
                    #launch pileup
                    snp.sam_pileup(ex,job,dictRun["url"],genomeRef,via=opt.via,stdout=pileupFilename)
                    
                    sampleName=job.groups[idGroup]['name']
                    if(nbRuns>1):
                        sampleName+=group['run_names'].get(idRun,str(idRun))
                    
                    dictPileupFile[pileupFilename]=sampleName

                #formate pileup file
                if(nbRuns>1):
                    # stat and correlation
                    print "many runs, need statistics"
 
            posAllUniqSNPFile=snp.posAllUniqSNP(ex,job,dictPileupFile)
            ex.add(posAllUniqSNPFile)

            for k in dictPileupFile.keys():
                ex.add(k)

            formatedPileupFilename=snp.parse_pileupFile(ex,job,dictPileupFile,posAllUniqSNPFile,via=opt.via)
            description="SNP analysis for samples: "+", ".join(dictPileupFile.values())
            description=set_file_descr("allSNP.txt",step="SNPs",type="txt")
            ex.add(formatedPileupFilename,description=description)
        
        allfiles = common.get_files(ex.id, M)
        logfile.close()
        print json.dumps(allfiles)

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())


