#!/bin/env python
"""
pickle_mappings.py
by Julien Delafontaine, <julien.delafontaine@epfl.ch>
09.08.2011

Takes a species identifier for Ensembl (e.g. 'human') and returns a pickle containing a tuple
(gene_ids, gene_names, trans_to_gene_mapping, exons_to_trans_mapping).
gene_ids is a list of gene IDs; gene_names is a dictionary {gene ID: gene name};
trans_to_gene_mapping is a dictionary {gene ID: [transcripts IDs]};
exons_to_trans_mapping is a dictionary {transcript ID: [exons IDs]}.
"""
import cogent.db.ensembl as ensembl
from numpy import *
import pickle, os, sys, urllib, time, random, math, numpy, getopt
from numpy.linalg import pinv

usage = """Usage: pickle_mappings.py [-h] species [-o output.pickle] [-r: ensembl_release]
>>> pickle_mappings.py human maptest -r 63

species: 
-h     Print this message and exit
-o     Name of the output file (default: species name)
-r     Version number of Ensembl release (default: 63)
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def generate_mappings(species='human',ensembl_release=63):
    ''' Return a tuple
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_by_gene, exons_by_trans).
    gene_ids is a list of gene IDs;
    gene_names is a dictionary         {gene ID: gene name};
    transcript_mapping is a dictionary {transcript ID: gene ID};
    exon_mapping is a dictionary       {exon ID: (transcript ID, gene ID)};
    trans_by_gene is a dictionary      {gene ID: [transcripts IDs]};
    exons_by_trans is a dictionary     {transcript ID: [exons IDs]}.

    - species is an Ensembl identifier as "human";
    - ensembl_release is the version number of the Ensembl release
    to extract data from.
    '''
    genome = ensembl.Genome(species, Release=ensembl_release)
    if genome.Species == 'None':
        print ensembl.Species
        raise Usage("Error: Unknown species")
    #genome.getDistinct('BioType')
    gene_ids = []; gene_names = {}
    transcript_mapping = {}; exon_mapping = {}
    trans_by_gene = {}; exons_by_trans = {}
    gene_generator = genome.getGenesMatching(BioType='protein_coding')
    #k = 0
    for g in gene_generator:
        #k+=1
        #if k>10: break;
        gid = g.StableId
        transcripts = g.Transcripts
        gene_ids.append(gid)
        gene_names[gid] = g.Symbol
        trans_by_gene[gid] = []
        for t in transcripts:
            tid = t.StableId
            exons = t.Exons
            trans_by_gene[gid].append(tid)
            transcript_mapping[tid] = gid
            exons_by_trans[tid] = []
            for e in exons:
                eid = e.StableId
                exons_by_trans[tid].append(eid)
                exon_mapping[eid] = (tid, gid)
    return (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_by_gene, exons_by_trans)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    ensembl_release = 63

    try:
        try:
            opts, args = getopt.getopt(argv, "h", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        
        if len(args) < 1:
            raise Usage("pickle_exon_mapping.py takes at least 1 argument.")

        species = str(args[0])
        output_filename = species
        
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            if o in ("-r",):
                ensembl_release = a
            if o in ("-o",):
                output_filename = a
            else:
                raise Usage("Unhandled option: " + o)

        ### Program body ###
        t1 = time.time()
        mappings = generate_mappings(species,ensembl_release)
        t2 = time.time()
        running_time = t2-t1
        print "Execution time:", str(running_time), "s."
        with open(output_filename+".pickle","wb") as f:
            pickle.dump(mappings,f)
        t3 = time.time()
        running_time = t3-t2
        print "Dumping time:", str(running_time), "s."

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    

if __name__ == '__main__':
    sys.exit(main())
