#!/bin/env python
"""
Interface to genrep functionalities.

"""

from bbcflib import genrep
from bbcflib.common import normalize_url
import sys, os, re
import optparse

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main():
    try:
        # Parse args
        usage = "genrep4humans.py [-h] [-l] [-f] [-b] -a [-t] [-r] [-u]"
        parser = optparse.OptionParser(usage=usage, description="Interface to genrep functionalities.")
        
        parser.add_option("-l","--list", action="store_true", default=False,
                          help="list available assemblies")
        parser.add_option("-f","--fasta", action="store_true", default=False,
                          help="get path to fasta")
        parser.add_option("-b","--bowtie", action="store_true", default=False,
                          help="get path prefix to bowtie indexes")
        parser.add_option("-s","--stats", action="store_true", default=False,
                          help="genome stats")
        parser.add_option("-i","--intype", type="int", default=0,
                          help="0 for genome (default), 1 for exonome, 2 for transcriptome")
        parser.add_option("-a","--assembly", default=None,
                          help="assembly (name or id)")
        parser.add_option("-t","--regions", default=None,
                          help="extract region to fasta (ex: 'chr2:2-1356;chr1:10-45')")
        parser.add_option("-o","--output", default=None,
                          help="output file (default standard output)")
        parser.add_option("-r","--root", default='/db/genrep/',
                          help="genrep root directory (default: '/db/genrep/')")
        parser.add_option("-u","--url", default='http://bbcftools.vital-it.ch/genrep/',
                          help="url to genrep (default: 'http://bbcftools.vital-it.ch/genrep/')")

        # Get variables
        (opt, args) = parser.parse_args()
        assembly_id = re.search('([._\-\w]+)', str(opt.assembly)).groups()[0]
        genrep_root = os.path.abspath(opt.root)
        genrep_url = normalize_url(opt.url)
        if opt.output:
            fout = open(re.search('([._\-\w]+)', str(opt.output)).groups()[0], 'w')
        else: fout = sys.stdout
        if opt.regions:
            regions = []; start=None; end=None
            for x in str(opt.regions).split(";"):
                chr,start,end = re.search('(\S+):(\d+)\-(\d+)',x).groups()[0:3]
                regions.append([chr,int(start),int(end)])

        # Program body
        g_rep = genrep.GenRep(url=genrep_url, root=genrep_root, intype=opt.intype)
        if opt.list:
            fout.write("\n".join(g_rep.assemblies_available())+"\n")
        if opt.assembly:
            assembly = g_rep.assembly(assembly_id)
        else:
            sys.exit(0)
        if opt.regions:
            seq = g_rep.fasta_from_regions(assembly.chromosomes, regions=regions, out={})[0]
            for reg in regions:
                write(">"+assembly.name+"|"+reg[0]+":"+str(reg[1])+"-"+str(reg[2])+"\n")
                write(seq[reg[0]].pop(0)+"\n")
        if opt.bowtie:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" bowtie index prefix\n")
            fout.write(assembly.index_path+"\n")
        if opt.fasta:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" fasta file\n")
            fout.write(g_rep.fasta_path(assembly)+"\n")
        if opt.stats:
            stats = g_rep.statistics(assembly,frequency=True)
            fout.write("\n".join([k+":\t"+str(stats[k]) for k in sorted(stats.keys())])+"\n")
        fout.close()
        
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    
if __name__ == '__main__':
    sys.exit(main())
