#!/usr/bin/env python

"""
Count reads on genes and transcripts from a genome-level BAM file and a
GTF file describing the exons, such as those provided by Emsembl or GenRep.
The GTF is assumed to be sorted at least w.r.t. chromosome name,
and the chromosome identifiers in the GTF must be the same as the BAM references.

Usage:
   rnacounter  [-n <int>] [-l <int>] [-s] [-m] [-c <string>] [-o <string>] BAM GTF
               [--version] [-h]

Options:
   -n <int>, --normalize <int>          Normalization constant [default: total number of reads].
   -l <int>, --fraglength <int>         Average fragment length [default: 350].
   -s, --stranded                       Compute sense and antisense reads separately [default: False].
   -m, --multiple                       Divide count by NH flag for multiply mapping reads [default: False].
   -c <string>, --chromosomes <string>  Chromosome names (comma-separated list).
   -o <string>, --output <string>       Output file to redirect stdout.
   -v, --version                        Displays version information and exits.
   -h, --help                           Displays usage information and exits.
"""

from rnacounter import rnacounter_main

from docopt import docopt
import os,sys

if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    bamname = os.path.abspath(args['BAM'])
    annotname = os.path.abspath(args['GTF'])
    if args['--output'] is None: output = sys.stdout
    if args['--chromosomes'] is None: chromosomes = []
    else: chromosomes = args['--chromosomes'].split(',')

    rnacounter_main(bamname,annotname,
                  multiple=args['--multiple'], stranded=args['--stranded'],
                  output=args['--output'], normalize=args['--normalize'],
                  chromosomes=chromosomes, fraglength=args['--fraglength'])

