#!/usr/bin/env python

"""
SNP detection workflow.
"""
import os,sys
from bbcflib import snp, mapseq, genrep
from run import run

class SNP_job:
    def __init__(self):
        self.opts = (("-v", "--via", "Run executions locally or using bsub (can be either 'local' or 'lsf')",
                                              {'default': "lsf"}),
                     ("-k", "--key", "Alphanumeric key of the new SNP job", {'default': None}),
                     ("-m", "--mapseq_limspath", "MiniLIMS where a previous Mapseq execution and files has been stored.",
                                              {'default': "/srv/mapseq/public/data/mapseq_minilims"}),
                     ("-w", "--working-directory", "Create execution working directories in wdir",
                                              {'default': os.getcwd(), 'dest':"wdir"}),
                     ("-c", "--config", "Config file", {'default': None}),
                     ("-d", "--minilims", "MiniLIMS where SNP executions and files will be stored.", \
                                              {'default': "/srv/snp/public/data/snp_minilims"}),
                     ("-f", "--fasta_path", "Path to a directory containing a fasta file for each chromosome",
                                              {'default':None}),
                     ("--mincov", "Minimum number of reads supporting an SNP at a position for it to be considered.",
                                              {'default':5}),
                     ("--minsnp", "Minimum percentage of reads supporting the SNP for it to be returned.",
                                              {'default':40}),
                    )
        self.name = "SNP job"
        self.usage = "run_snp.py [OPTIONS]"
        self.desc = """Compares sequencing data to a reference assembly to detect SNP."""
        self.hts = 'hts_snp'

    def workflow(self,ex,job,opt,gl,logfile,debugfile):
        mapseq_url = gl.get('hts_mapseq',{}).get('url')
        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )
        logfile.write("Fetch bam and wig files.\n");logfile.flush()
        job = mapseq.get_bam_wig_files(ex, job, minilims=opt.mapseq_limspath, hts_url=mapseq_url,
                                       script_path=gl.get('script_path',''), via=opt.via)
        mincov = job.options.get('mincov') or opt.mincov
        minsnp = job.options.get('minsnp') or opt.minsnp
        snp.snp_workflow(ex,job,assembly,mincov=mincov,minsnp=minsnp,path_to_ref=opt.fasta_path,via=opt.via)

if __name__ == '__main__':
    sys.exit(run(SNP_job()))



