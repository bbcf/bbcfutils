#!/usr/bin/env python

"""
Command-line interface to genrep functionalities.

"""
from bbcflib import genrep
from bbcflib.common import normalize_url
import sys, os, re
import optparse

usage = "genrep4humans.py [OPTIONS]"
descr = "Command-line interface to genrep functionalities."
opts = (("-l", "--list", "list available assemblies, or a chromosome table if an assembly is specified",
         {'action': "store_true", 'default': False}),
        ("-f", "--fasta", "get path to fasta", {'action': "store_true", 'default': False}),
        ("-b", "--bowtie", "get path prefix to bowtie indexes", {'action': "store_true", 'default': False}),
        ("-B", "--bowtie2", "get path prefix to bowtie2 indexes", {'action': "store_true", 'default': False}),
        ("-d", "--db", "get path to sqlite database", {'action': "store_true", 'default': False}),
        ("-s", "--stats", "genome stats", {'action': "store_true", 'default': False}),
        ("-i", "--intype", "0 for genome (default), 1 for exonome, 2 for transcriptome", {'type': "int", 'default': 0}),
        ("-a", "--assembly", "assembly (name or id)", {'default': None}),
        ("-t", "--regions", "extract regions to fasta (a string: 'chr2:2-1356,chr1:10-45' or a bed/gff/sql filename)", {'default': None}),
        ("-o", "--output", "output file (default: standard output)", {'default': None}),
        ("-r", "--root", "genrep root directory (default: '/db/genrep/')", {'default': '/db/genrep/'}),
        ("-u", "--url", "url to genrep (default: 'http://bbcf-serv01.epfl.ch/genrep/')",{'default': 'http://bbcf-serv01.epfl.ch/genrep/'}),
        ("-g", "--genes", "extract coordinates for a list of Ensembl ids, specified as a comma-separated list or a filename (gene id must be the first word of each line in the file)",{}),
        ("-z", "--all", "extract coordinates of all genes/transcripts/exons (depending on the 'intype')",{'action': "store_true", 'default': False}),
        ("-c", "--convert", "convert bam headers to natural chromosome names",{}))

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def _parse_list(filename):
    with open(filename) as fin:
        for row in fin:
            gid = re.search(r'\S+',row)
            if gid:
                yield gid.group()

def main():
    try:
        # Parse args
        parser = optparse.OptionParser(usage=usage, description=descr)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])

        # Get variables
        (opt, args) = parser.parse_args()
        if opt.assembly:
            assembly_id = re.search('([._\-\w]+)', str(opt.assembly)).groups()[0]
        genrep_root = os.path.abspath(opt.root)
        genrep_url = normalize_url(opt.url)
        if opt.output:
            fout = open(re.search('([._\-\w]+)', str(opt.output)).groups()[0], 'w')
        else:
            fout = sys.stdout
        regions = None
        if opt.regions:
            if os.path.exists(opt.regions):
                regions = opt.regions
            else:
                regions = []
                for x in str(opt.regions).split(","):
                    chrom,start,end = re.search('(\S+):(\d+)\-(\d+)',x).groups()[0:3]
                    regions.append([chrom,int(start),int(end)])

        # Program body
        g_rep = genrep.GenRep(url=genrep_url, root=genrep_root)
        if opt.assembly:
            assembly = genrep.Assembly(assembly=assembly_id,genrep=g_rep,intype=opt.intype)
        if opt.list:
            if opt.assembly:
                table = ["\t".join((v['ac'],k,str(v['length'])))
                         for k,v in assembly.chrmeta.iteritems()]
                fout.write("\n".join(table)+"\n")
            else:
                fout.write("\n".join(v[1] for v in g_rep.assemblies_available())+"\n")
            return 0
        if not(opt.assembly):
            parser.print_help()
            return 0
        if regions:
            seq = assembly.fasta_from_regions(regions=regions, out=fout)[0]
        if opt.bowtie:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" bowtie index prefix\n")
            fout.write(assembly.index_path+"\n")
        if opt.bowtie2:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" bowtie2 index prefix\n")
            fout.write(re.sub(r'bowtie/','bowtie2/',assembly.index_path)+"\n")
        if opt.fasta:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" fasta file\n")
            fout.write(assembly.fasta_path()+"\n")
        if opt.db:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" sqlite file\n")
            fout.write(assembly.sqlite_path+"\n")
        if opt.genes:
            if os.path.exists(opt.genes):
                glist = _parse_list(opt.genes)
            else:
                glist = opt.genes.split(",")
            for gcoord in assembly.gene_coordinates(glist):
                fout.write("\t".join([str(x) for x in gcoord])+"\n")
        if opt.all:
            if opt.intype == 1:
                feats = assembly.exon_track()
            elif opt.intype == 2:
                feats = assembly.transcript_track()
            else:
                feats = assembly.gene_track()
            for gcoord in feats:
                fout.write("\t".join([str(x) for x in gcoord])+"\n")
        if opt.stats:
            stats = assembly.statistics(frequency=True)
            fout.write("\n".join([k+"\t"+str(stats[k]) for k in sorted(stats.keys())])+"\n")
        fout.close()
        if opt.convert:
            if not(os.path.exists(opt.convert)):
                raise Usage("No such file: %s."%opt.convert)
            if not(opt.output):
                raise Usage("Need an output file name.")
            import pysam
            infile = pysam.Samfile( opt.convert )
            header = infile.header
            chromosomes = dict((v['ac'],k) for k,v in assembly.chrmeta.iteritems())
            for h in header["SQ"]:
                if h["SN"] in chromosomes:
                    h["SN"] = chromosomes[h["SN"]]
            outfile = pysam.Samfile(re.search('([._\-\w]+)', str(opt.output)).groups()[0], 'wb', header=header )
            for read in infile:
                outfile.write(read)
            outfile.close()
            infile.close()

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
