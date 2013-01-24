#!/usr/bin/env python

"""
A High-throughput RNA-seq analysis workflow.

python run_rnaseq.py -v lsf -c config_files/gapkowt.txt -d rnaseq -p transcripts,genes
python run_rnaseq.py -v lsf -c config_files/rnaseq.txt -d rnaseq -p genes -m /scratch/cluster/monthly/jdelafon/mapseq
"""
import os,sys
from bbcflib import rnaseq, mapseq
from run import run

class Rnaseq_job:
    def __init__(self):
        self.opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default':"lsf"}),
                     ("-k", "--key", "Alphanumeric key of the new RNA-seq job", {'default': None}),
                     ("-d", "--minilims", "MiniLIMS where RNAseq executions and files will be stored.",
                                              {'default': "/srv/rnaseq/public/data/rnaseq_minilims"}),
                     ("-m", "--mapseq_minilims", "MiniLIMS where a previous Mapseq execution and files has been stored.",
                                              {'default': "/srv/mapseq/public/data/mapseq_minilims"}),
                     ("-w", "--working-directory", "Create execution working directories in wdir",
                                              {'default': os.getcwd(), 'dest':"wdir"}),
                     ("-c", "--config", "Config file", {'default': None}),
                     ("-p", "--pileup_level", "Target features, inside of quotes, separated by commas.\
                                              E.g. 'genes,exons,transcripts'",{'default':"genes,exons,transcripts"}),
                     ("-j", "--junctions", "Whether or not to seqrch for splice junctions using SOAPsplice",
                                              {'action':"store_true", 'default':False}),
                    )
        self.name = "RNA-seq job"
        self.usage = "run_rnaseq.py [-h -v via -k key -c config_file -w working_directory -d minilims -m mapseq_minilims]"
        self.desc = """A High-throughput RNA-seq analysis workflow. It returns a file containing
                  a column of transcript counts for each given BAM file, normalized using DESeq's
                  size factors. """
        self.hts = 'hts_rnaseq'

    def workflow(self,ex,job,opt,gl,logfile,debugfile):
        job.options['discard_pcr_duplicates'] = False
        job.options['find_junctions'] = opt.junctions or (job.options.get('find_junctions','').lower() in ['1','true','t'])
        mapseq_url = gl.get('hts_mapseq',{}).get('url')
        pileup_level = opt.pileup_level.split(',')
        logfile.write("Fetch bam and wig files\n");logfile.flush()
        job = mapseq.get_bam_wig_files(ex, job, minilims=opt.mapseq_minilims, hts_url=mapseq_url,
                                       script_path=gl.get('script_path',''), via=opt.via, fetch_unmapped=True)
        logfile.write("Starting workflow.\n");logfile.flush()
        rnaseq.rnaseq_workflow(ex, job, pileup_level=pileup_level, via=opt.via,
                               rpath=gl.get('script_path'), junctions=job.options['find_junctions'],
                               debugfile=debugfile, logfile=logfile)

if __name__ == '__main__':
    sys.exit(run(Rnaseq_job()))


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#
