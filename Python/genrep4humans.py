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
        ("-s", "--stats", "genome stats", {'action': "store_true", 'default': False}),
        ("-i", "--intype", "0 for genome (default), 1 for exonome, 2 for transcriptome", {'type': "int", 'default': 0}),
        ("-a", "--assembly", "assembly (name or id)", {'default': None}),
        ("-t", "--regions", "extract regions to fasta (ex: 'chr2:2-1356,chr1:10-45')", {'default': None}),
        ("-o", "--output", "output file (default standard output)", {'default': None}),
        ("-r", "--root", "genrep root directory (default: '/db/genrep/')", {'default': '/db/genrep/'}),
        ("-u", "--url", "url to genrep (default: 'http://bbcftools.vital-it.ch/genrep/')",{'default': 'http://bbcftools.vital-it.ch/genrep/'}),
        ("-c", "--convert", "convert bam headers to natural chromosome names",{}))

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def _compact_key(key):
    if key == None or isinstance(key,str):
        return key
    elif isinstance(key,tuple) and len(key)>2:
        return str(key[0])+"_"+str(key[1])+"."+str(key[2])
    else:
        raise ValueError("Can't handle this chromosomes key ",key)

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
        if opt.regions:
            regions = []; start=None; end=None
            for x in str(opt.regions).split(","):
                chrom,start,end = re.search('(\S+):(\d+)\-(\d+)',x).groups()[0:3]
                regions.append([chrom,int(start),int(end)])

        # Program body
        g_rep = genrep.GenRep(url=genrep_url, root=genrep_root)
        if opt.assembly:
            assembly = genrep.Assembly(assembly=assembly_id,genrep=g_rep,intype=opt.intype)
        if opt.list:
            if opt.assembly:
                table = ["\t".join((_compact_key(k),v['name'],str(v['length'])))
                         for k,v in assembly.chromosomes.iteritems()]
                fout.write("\n".join(table)+"\n")
            else:
                fout.write("\n".join(g_rep.assemblies_available())+"\n")
            return 0
        if not(opt.assembly):
            return 0
        if opt.regions:
            seq = assembly.fasta_from_regions(regions=regions, out={})[0]
            for reg in regions:
                fout.write(">"+assembly.name+"|"+reg[0]+":"+str(reg[1])+"-"+str(reg[2])+"\n")
                fout.write(seq[reg[0]].pop(0)+"\n")
        if opt.bowtie:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" bowtie index prefix\n")
            fout.write(assembly.index_path+"\n")
        if opt.fasta:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" fasta file\n")
            fout.write(assembly.fasta_path()+"\n")
        if opt.stats:
            stats = assembly.statistics(frequency=True)
            fout.write("\n".join([k+":\t"+str(stats[k]) for k in sorted(stats.keys())])+"\n")
        fout.close()
        if opt.convert:
            if not(os.path.exists(opt.convert)):
                raise ValueError("No such file: %s."%opt.convert)
            if not(opt.output):
                raise ValueError("Need an output file name.")
            import pysam
            infile = pysam.Samfile( opt.convert, 'rb' )
            header = infile.header
            chromosomes = dict((_compact_key(k),v['name']) for k,v in assembly.chromosomes.iteritems())
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
