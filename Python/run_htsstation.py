#!/usr/bin/env python
"""
Top-level script for HTSstation.epfl.ch pipelines.

"""

import optparse
from bbcflib.workflows import Workflow, Usage

_module_list = ["demultiplexing","mapseq","chipseq","rnaseq","snp","4cseq"]
_mapseq_lims_opt = ("-m", "--mapseq_minilims", 
                    "MiniLIMS where a previous Mapseq execution and files have been stored.",
                    {'default': None})

############################# Demultiplexing #############################
class DemulitplexWorkflow(Workflow):

    def __init__(self):
        opts = (,)
        usage = ""
        desc = "Demultiplexing routine for 4C-seq data."
        Workflow.__init__(self,module="demultiplex",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        Workflow.check_options(self)
        self.main_args = {"job": self.job,"gl": self.globals}
        return True

    def init_files(self,ex):
        from bbcflib.mapseq import get_fastq_files
        self.log_write("fetch fastq files.")
        self.job = get_fastq_files( ex, self.job, self.job.dafl )
        return None

############################### MapSeq ###############################
class MapseqWorkflow(Workflow):

    def __init__(self):
        opts = (("--noqc","Skip fastqc step",{'action':"store_true"}),)
        usage = "[--noqc]"
        desc = "A High-throughput sequencing data mapping workflow."
        Workflow.__init__(self,module="mapseq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        Workflow.check_options(self)
        map_args = self.job.options.get('map_args',{})
        self.job.options['create_gdv_project'] &= self.job.options['compute_densities']
        self.main_args = {"job": self.job,
                          "assembly": self.job.assembly,
                          "map_args": map_args,
                          "gl": self.globals,
                          "via": self.opt.via,,
                          "debugfile": self.debugfile, 
                          "logfile": self.logfile}

        return True

    def init_files(self,ex):
        self.log_write("fetch fastq files.")
        self.job = mapseq.get_fastq_files( ex, self.job, self.job.dafl )
        if not self.opt.noqc: 
            self.log_write("Generate QC report.")
            mapseq.run_fastqc( ex, job, via=self.opt.via )
        return None

############################### ChipSeq ###############################
class ChipseqWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,)
        usage = "[-m mapseq_minilims]"
        desc = """A High-throughput ChIP-seq peak analysis workflow."""
        Workflow.__init__(self,module="chipseq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        Workflow.check_option(self)
        self.main_args = {"job": self.job,
                          "assembly": self.job.assembly,
                          "script_path": self.globals.get('script_path',''),
                          "logfile": self.logfile, 
                          "via": self.opt.via}
        return True

############################### RnaSeq ###############################
class RnaseqWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,
                ("-p", "--pileup_level", 
                 "Target features, inside of quotes, separated by commas. E.g. 'genes,exons,transcripts'",
                 {'default':"genes,exons,transcripts"}),
                ("-j", "--junctions", "whether or not to search for splice junctions using soapsplice",
                 {'action':"store_true", 'default':False}),
                ("-u", "--unmapped", "whether or not to use unmapped reads from a previous mapping job, if available",
                 {'action':"store_true", 'default':False}))
        usage = "[-m mapseq_minilims] [-p -j -u]"
        desc = """A High-throughput RNA-seq analysis workflow. It returns text files containing
                  read counts for exons, genes and transcripts for each given BAM file, and the result
                  of a DESeq run (differential expression analysis) for every pair of groups. """
        Workflow.__init__(self,module="rnaseq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        def = {'discard_pcr_duplicates': (False,),
               'find_junctions': (False,),
               'unmapped': (False,)}
        Workflow.check_option(self, def)
        self.main_args = {"job": self.job, 
                          "pileup_level": opt.pileup_level.split(','),
                          "via": self.opt.via,
                          "rpath": self.globals.get('script_path',''), 
                          "junctions": job.options['find_junctions'], 
                          "unmapped": job.options['unmapped'],
                          "debugfile": self.debugfile, 
                          "logfile": self.logfile}
        return True

############################### Snps ###############################
class SnpWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,
                ("-f", "--fasta_path", "Path to a directory containing a fasta file for each chromosome",
                 {'default':None}),
                ("--mincov", "Minimum number of reads supporting an SNP at a position for it to be considered.",
                 {'default':5}),
                ("--minsnp", "Minimum percentage of reads supporting the SNP for it to be returned.",
                 {'default':40}))
        usage = "[-m mapseq_minilims] [--mincov --minsnp]"
        desc = """Compares sequencing data to a reference assembly to detect SNPs."""
        Workflow.__init__(self,module="snp",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        def = {}
        Workflow.check_option(self, def)
        mincov = job.options.get('mincov') or opt.mincov
        minsnp = job.options.get('minsnp') or opt.minsnp
        self.main_args = {"job": self.job, 
                          "assembly": self.job.assembly,
                          "mincov": mincov,
                          "minsnp": minsnp,
                          "path_to_ref": opt.fasta_path,
                          "via": self.opt.via}
        return True

############################### 4C-seq ###############################
class C4seqWorkflow(Workflow):

    def __init__(self):
        opts = (_mapseq_lims_opt,)
        usage = "[-m mapseq_minilims]"
        desc = """A High-throughput 4C-seq analysis workflow."""
        Workflow.__init__(self,module="c4seq",opts=opts,usage=usage,desc=desc)

    def check_options(self):
        Workflow.check_option(self)
        primers_dict = c4seq.loadPrimers(os.path.join(self.opt.wdir,'primers.fa'))
        self.main_args = {"job": self.job, 
                          "primers_dict": primers_dict,
                          "assembly": self.job.assembly,
                          "c4_url": self.globals.get('hts_4cseq',{}).get('url'),
                          "script_path": self.globals.get('script_path',''),
                          "logfile": self.logfile, 
                          "via": self.opt.via}
        return True

############################### MAIN ###############################
def main():
    try:
        module = None
        if len(sys.argv) > 1: module = sys.argv.pop(1)
        if module not in _module_list:
            m = module and "No such operation: %s, choose one of %s." %(module,str(_module_list)) or ''
            raise Usage(m)
        if module == "demultiplexing":
            WF = DemulitplexWorkflow()
        elif module == "mapseq":
            WF = MapseqWorkflow()
        elif module == "chipseq":
            WF = ChipseqWorkflow()
        elif module == "rnaseq":
            WF = RnaseqWorkflow()
        elif module == "snp":
            WF = SnpWorkflow()
        elif module == "4cseq":
            WF = C4seqWorkflow()

        parser = optparse.OptionParser(usage=WF.usage, description=WF.desc)
        for opt in WF.opts:
            if len(opt) == 4:
                parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
            elif len(opt) == 3:
                parser.add_option(opt[0],help=opt[1],**opt[2])
        (opt, args) = parser.parse_args()
        return WF(opt)

    except Usage, err:
        print >>sys.stderr, err.msg
        if parser: parser.print_help()
        return 2

if __name__ == '__main__':
    sys.exit(main())
