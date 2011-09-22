#!/bin/env python
"""
Interface to genrep functionalities.

"""
from bbcflib import genrep
from bbcflib.common import normalize_url
import sys, getopt, os, re

usage = """genrep4humans.py [-h] [-l] [-f] [-b] -a [-t] [-r] [-u]
-h --help     print this message and exit
-l --list     list available assemblies
-f --fasta    get path to fasta
-b --bowtie   get path prefix to bowtie indexes
-s --stats    genome stats
-a --assembly assembly (name or id)
-r --regions  extract region to fasta (ex: 'chr2:2-1356;chr1:10-45')
-t --type     data type (0=genome, 1=exome, 2=transcriptome, default 0)
-o --output   output file (default standard output)
-r --root     genrep root directory (default: '/db/genrep/')
-u --url      url to genrep (default: 'http://bbcftools.vital-it.ch/genrep/')
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    genrep_url = 'http://bbcftools.vital-it.ch/genrep/'
    genrep_root = '/db/genrep/'
    assembly_id = -1
    intype = 0
    fasta = False
    bowtie = False
    stats = False
    list = False
    fout = sys.stdout
    regions = None
    start = None
    end = None
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hlfbsa:r:t:o:r:u:",
                                      ["help","list","fasta","bowtie","stats",
                                       "type=","assembly=","regions=","output=",
                                       "root=","url="])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                return 0
            elif o in ("-l", "--list"):
                list = True
            elif o in ("-f", "--fasta"):
                fasta = True
            elif o in ("-b", "--bowtie"):
                bowtie = True
            elif o in ("-s", "--stats"):
                stats = True
            elif o in ("-a", "--assembly"):
                assembly_id = re.search('([._\-\w]+)',a).groups()[0]
            elif o in ("-r", "--regions"):
                regions = []
                for x in a.split(";"):
                    chr,start,end = re.search('(\S+):(\d+)\-(\d+)',x).groups()[0:3]
                    regions.append([chr,int(start),int(end)])
            elif o in ("-t", "--type"):
                intype = int(a)
            elif o in ("-o", "--output"):
                fout = open(re.search('([._\-\w]+)',a).groups()[0], 'w')
            elif o in ("-r", "--root"):
                genrep_root = os.path.abspath(a)
            elif o in ("-u", "--url"):
                genrep_url = normalize_url(a)
            else:
                raise Usage("Unhandled option: " + o)
        g_rep = genrep.GenRep( url=genrep_url, root=genrep_root, intype=intype )
        if list:
            fout.write("\n".join(g_rep.assemblies_available())+"\n")
        if assembly_id>0:
            assembly = g_rep.assembly( assembly_id )
        else:
            return 0
        if regions:
            seq = g_rep.fasta_from_regions(assembly.chromosomes,regions=regions,out={})[0]
            for reg in regions:
                fout.write(">"+assembly.name+"|"+reg[0]+":"+str(reg[1])+"-"+str(reg[2])+"\n")
                fout.write(seq[reg[0]].pop(0)+"\n")
        if bowtie:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" bowtie index prefix\n")
            fout.write(assembly.index_path+"\n")
        if fasta:
            fout.write(">"+str(assembly.id)+":"+assembly.name+" fasta file\n")
            fout.write(g_rep.fasta_path(assembly)+"\n")
        if stats:
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
