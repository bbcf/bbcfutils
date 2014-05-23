#!/usr/bin/env python
"""
Top-level script for HTSstation.epfl.ch pipelines.

"""

import optparse, sys, os
from bbcflib.workflows import Workflow, Usage

_module_list = ["demultiplexing","mapseq","chipseq","rnaseq","snp","4cseq"]
_mapseq_lims_opt = ("-m", "--mapseq_minilims",
                    "MiniLIMS where a previous Mapseq execution and files have been stored.",
                    {'default': None})

############################# Demultiplexing #############################
class DemulitplexWorkflow(Workflow):

    def __init__(self):
        desc = "Demultiplexing routine for 4C-seq data."
        Workflow.__init__(self,module="demultiplex",name="demultiplexing",desc=desc)

    def check_options(self):
        Workflow.check_options(self)
        self.main_args = {"job": self.job,
                          "gl": self.globals,
                          "via": self.opts.via,
                          "debugfile": self.debugfile,
                          "logfile": self.logfile}
        return True

    def init_files(self,ex):
        from bbcflib.mapseq import get_fastq_files
        self.log_write("fetch fastq files.")
        self.job = get_fastq_files( ex, self.job )
        return None

############################### MapSeq ###############################
class MapseqWorkflow(Workflow):

    def __init__(self):
        opts = (("--noqc","Skip fastqc step",{'action':"store_true"}),)
        usage = "[--noqc]"
        desc = "A High-throughput sequencing data mapping workflow."
        Workflow.__init__(self,module="mapseq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        Workflow.check_options(self,{'bowtie2': (True,)})
        map_args = self.job.options.get('map_args',{})
        self.job.options['create_gdv_project'] &= self.job.options['compute_densities']
        self.main_args = {"job": self.job,
                          "assembly": self.job.assembly,
                          "map_args": map_args,
                          "gl": self.globals,
                          "bowtie2": self.job.options['bowtie2'],
                          "via": self.opts.via,
                          "debugfile": self.debugfile,
                          "logfile": self.logfile}

        return True

    def init_files(self,ex):
        self.log_write("fetch fastq files.")
        self.job = self.sysmod.get_fastq_files( ex, self.job )
        if not self.opts.noqc:
            self.log_write("Generate QC report.")
            self.sysmod.run_fastqc( ex, self.job, via=self.opts.via )
        return None

############################### ChipSeq ###############################
class ChipseqWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,)
        usage = "[-m mapseq_minilims]"
        desc = """A High-throughput ChIP-seq peak analysis workflow."""
        Workflow.__init__(self,module="chipseq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        Workflow.check_options(self)
        self.main_args = {"job_or_dict": self.job,
                          "assembly": self.job.assembly,
                          "script_path": self.globals.get('script_path',''),
                          "logfile": self.logfile,
                          "via": self.opts.via}
        return True

############################### RnaSeq ###############################
class RnaseqWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,
                ("-p", "--pileup_level",
                 "Target features, inside of quotes, separated by commas. E.g. 'genes,exons,transcripts'.",
                 {'default':"genes,exons,transcripts"}),
                ("-j", "--junctions", "Whether or not to search for splice junctions using soapsplice.",
                 {'action':"store_true", 'default':False}),
                ("--stranded", "To indicate that the protocol was strand-specific.",
                 {'action':"store_true", 'default':False}),
               )
        usage = "[-m mapseq_minilims] [-p -j --stranded]"
        desc = """A High-throughput RNA-seq analysis workflow. It returns text files containing
                  read counts for exons, genes and transcripts for each given BAM file, and the result
                  of a DESeq run (differential expression analysis) for every pair of groups. """
        Workflow.__init__(self,module="rnaseq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        more_defs = {'discard_pcr_duplicates': (False,),
                     'find_junctions': (False,),
                     'stranded': (False,)}
        Workflow.check_options(self, more_defs)
        self.main_args = {"job": self.job,
                          "assembly": self.job.assembly,
                          "pileup_level": self.opts.pileup_level.split(','),
                          "via": self.opts.via,
                          "rpath": self.globals.get('script_path',''),
                          "juliapath": self.globals.get('script_path',''),
                          "junctions": self.job.options['find_junctions'] or self.opts.junctions,
                          "stranded": self.job.options['stranded'] and not self.opts.stranded,
                          "debugfile": self.debugfile,
                          "logfile": self.logfile}
        return True

############################### Snps ###############################
class SnpWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,
                ("-f", "--fasta_path", "Path to a directory containing a fasta file for each chromosome",
                 {'default':None}),
                ("--mincov", "Minimum coverage to call the SNP.",
                 {'default':5}),
                ("--minsnp", "Minimum percentage of reads per allele to call the SNP.",
                 {'default':40}))
        usage = "[-m mapseq_minilims] [--mincov --minsnp]"
        desc = """Compares sequencing data to a reference assembly to detect SNPs."""
        Workflow.__init__(self,module="snp",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        more_defs = {'bowtie2': (True,)}
        Workflow.check_options(self, more_defs)
        mincov = int(self.job.options.get('mincov') or self.opts.mincov)
        minsnp = int(self.job.options.get('minsnp') or self.opts.minsnp)
        self.main_args = {"job": self.job,
                          "assembly": self.job.assembly,
                          "mincov": mincov,
                          "minsnp": minsnp,
                          "path_to_ref": self.opts.fasta_path,
                          "via": self.opts.via}
        return True

############################### 4C-seq ###############################
class C4seqWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,)
        usage = "[-m mapseq_minilims]"
        desc = """A High-throughput 4C-seq analysis workflow."""
        Workflow.__init__(self,module="c4seq",name="4cseq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        Workflow.check_options(self)
        self.suffix = ['merged']
        primers_dict = self.sysmod.loadPrimers(os.path.join(self.opts.wdir,'primers.fa'))
        self.main_args = {"job": self.job,
                          "primers_dict": primers_dict,
                          "assembly": self.job.assembly,
                          "script_path": self.globals.get('script_path',''),
                          "logfile": self.logfile,
                          "via": self.opts.via}
        return True

############################### MAIN ###############################
def main():
    parser = None
    try:
        module = ''
        if len(sys.argv) > 1: module = sys.argv.pop(1)
        if module[:4] == "demu":
            WF = DemulitplexWorkflow()
        elif module[:3] == "map":
            WF = MapseqWorkflow()
        elif module[:4] == "chip":
            WF = ChipseqWorkflow()
        elif module[:3] == "rna":
            WF = RnaseqWorkflow()
        elif module[:3] == "snp":
            WF = SnpWorkflow()
        elif module[:2] in ["c4","4c"]:
            WF = C4seqWorkflow()
        else:
            m = "No such operation: %s, choose one of %s." %(module,str(_module_list)) if module else ''
            raise Usage(m)

        parser = optparse.OptionParser(usage=WF.usage, description=WF.desc)
        for opt in WF.opts:
            if len(opt) == 4:
                parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
            elif len(opt) == 3:
                parser.add_option(opt[0],help=opt[1],**opt[2])
        (opt, args) = parser.parse_args()
        return WF(opt)

    except Usage, err:
        print >>sys.stderr, '\n',err.msg,'\n'
        if parser: parser.print_help()
        else: print "run_htsstation.py %s [OPTIONS]" %str(_module_list)
        return 2

if __name__ == '__main__':
    sys.exit(main())
