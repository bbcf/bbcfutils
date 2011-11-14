#!/usr/bin/env python

"""
pickle_mappings.py
by Julien Delafontaine, <julien.delafontaine@epfl.ch>
09.08.2011

Takes a species identifier for Ensembl (e.g. 'human') and returns a pickle containing a tuple
(gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans).
gene_ids is a list of gene IDs;
gene_names is a dictionary         {gene ID: gene name};
transcript_mapping is a dictionary {transcript ID: gene ID};
exon_mapping is a dictionary       {exon ID: ([transcript IDs], gene ID)};
trans_in_gene is a dictionary      {gene ID: [transcripts IDs]};
exons_in_trans is a dictionary     {transcript ID: [exons IDs]}.
"""
import cogent.db.ensembl as ensembl
import pickle, sys, getopt

usage = """Usage: pickle_mappings.py [-h] [-o output.pickle] [-r ensembl_release] species
>>> pickle_mappings.py -o maptest -r 63 human

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
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans).
    gene_ids is a list of gene IDs;
    gene_names is a dictionary         {gene ID: gene name};
    transcript_mapping is a dictionary {transcript ID: gene ID};
    exon_mapping is a dictionary       {exon ID: ([transcript IDs], gene ID)};
    trans_in_gene is a dictionary      {gene ID: [transcripts IDs]};
    exons_in_trans is a dictionary     {transcript ID: [exons IDs]}.

    - species is an Ensembl identifier as "Homo Sapiens" or "C.Elegans";
    - ensembl_release is the version number of the Ensembl release
    to extract data from.
    '''
    genome = ensembl.Genome(species, Release=ensembl_release)
    if genome.Species == 'None':
        print >>sys.stderr, ensembl.Species
        raise Usage("Error: Unknown species")
    #genome.getDistinct('BioType')
    gene_ids = []; gene_names = {}
    transcript_mapping = {}; exon_mapping = {}
    trans_in_gene = {}; exons_in_trans = {}
    gene_generator = genome.getGenesMatching(BioType='protein_coding')
    #k = 0
    for g in gene_generator:
        #k+=1
        #if k>10: break;
        gid = str(g.StableId)
        gene_ids.append(gid)
        gene_names[gid] = str(g.Symbol)
        trans_in_gene[gid] = [str(t.StableId) for t in g.Transcripts]
        allexons = []
        for t in g.Transcripts:
            allexons.extend(t.Exons)
            transcript_mapping[str(t.StableId)] = gid
            exons_in_trans[str(t.StableId)] = [str(e.StableId) for e in t.Exons]
        for e in allexons:
            exon_mapping[str(e.StableId)] = ([str(t.StableId) for t in g.Transcripts
                                            if e.StableId in exons_in_trans[t.StableId]], gid)
    return (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    try:
        try:
            opts, args = getopt.getopt(argv, "ho:r:", ["help","output=","release="])
        except getopt.error, msg:
            raise Usage(msg)

        if len(args) < 1:
            raise Usage("pickle_exon_mapping.py takes at least 1 argument.")

        ensembl_release = 63
        species = str(args[0])
        output_filename = species

        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-r", "--release"):
                ensembl_release = a
            elif o in ("-o", "--output"):
                output_filename = a
            else:
                raise Usage("Unhandled option: " + o)

        ### Program body ###
        mappings = generate_mappings(species,ensembl_release)
        with open(output_filename+".pickle","wb") as f:
            pickle.dump(mappings,f)

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2


if __name__ == '__main__':
    sys.exit(main())
