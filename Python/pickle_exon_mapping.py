#!/bin/env python
"""
pickle_exon_mapping.py
by Fred Ross, <fred.ross@epfl.ch>
2011-05-16

Takes a FASTA file of exons and produces a pickle containing a 2-tuple (*gene labels*, *exon mapping*).  Each exon is assigned an arbitrary, unique, integer ID, as is each gene.  The *i*th entry of the list *gene labels* is a string giving the name of the gene with ID *i*.  The *i*th entry of *exon mapping* is an integer giving the ID of the gene that exon is part of.
"""
import os
import sys
import getopt
import pickle
from Bio import SeqIO

usage = """Usage: pickle_exon_mapping.py [-h] input.fast output.pickle

-h     Print this message and exit
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def generate_mapping(sequences):
    mapping = []
    label_to_id = {}
    id_to_label = {}
    next_gene_id = 0
    for i,s in enumerate(sequences):
        # label_to_id and id_to_label contain all genes used up to this point
        # mapping contains the gene id of the ith exon, up to this point
        gene_label = s.id.split('|')[1]
        if not(label_to_id.has_key(gene_label)):
            label_to_id[gene_label] = next_gene_id
            id_to_label[next_gene_id] = gene_label
            next_gene_id += 1
        mapping.append(label_to_id[gene_label])
    # After the loop, next_gene_id = # of genes
    gene_labels = [id_to_label[i] for i in xrange(next_gene_id)]
    return (gene_labels, mapping)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        try:
            opts, args = getopt.getopt(argv, "h", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            else:
                raise Usage("Unhandled option: " + o)

        if len(args) != 2:
            raise Usage("pickle_exon_mapping.py takes exactly two arguments.")

        elif not(os.path.exists(args[0])):
            raise Usage("Input file does not exist.")
        elif os.path.exists(args[1]):
            raise Usage("Output file already exists.")
        else:
            [input_filename, output_filename] = args

        # Program body
        seqs = SeqIO.parse(input_filename, 'fasta')
        d = generate_mapping(seqs)
        with open(output_filename, 'wb') as out:
            pickle.dump(d, out)

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    

if __name__ == '__main__':
    sys.exit(main())
