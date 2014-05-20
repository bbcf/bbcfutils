#!/usr/bin/env python

"""
Count reads on genes and transcripts from a genome-level BAM file and a
GTF file describing the exons, such as those provided by Emsembl or GenRep.
The GTF is assumed to be sorted at least w.r.t. chromosome name,
and the chromosome identifiers in the GTF must be the same as the BAM references.

Usage:
   rnacounter  [-t TYPE] [-n <int>] [-l <int>] [-s] [-m] [-c CHROMS] [-o OUTPUT] BAM GTF
               [--version] [-h]

Options:
   -t TYPE, --type TYPE                 Type of genomic feature to count on: 'genes' or 'transcripts'
                                        [default: genes].
   -n <int>, --normalize <int>          Normalization constant [default: total number of reads].
   -l <int>, --fraglength <int>         Average fragment length [default: 350].
   -s, --stranded                       Compute sense and antisense reads separately [default: False].
   -m, --multiple                       Divide count by NH flag for multiply mapping reads [default: False].
   -c CHROMS, --chromosomes CHROMS      Chromosome names (comma-separated list).
   -o OUTPUT, --output OUTPUT           Output file to redirect stdout.
   -v, --version                        Displays version information and exits.
   -h, --help                           Displays usage information and exits.
"""

from rnacounter import rnacounter_main
import os,sys


from docopt import docopt

if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    bamname = os.path.abspath(args['BAM'])
    annotname = os.path.abspath(args['GTF'])
    if args['--output'] is None: args['--output'] = sys.stdout
    else: args['--output'] = open(args['--output'], "wb")
    if args['--chromosomes'] is None: args['--chromosomes'] = []
    else: args['--chromosomes'] = args['--chromosomes'].split(',')
    assert args['--type'].lower() in ["genes","transcripts"], \
        "TYPE must be one of 'genes' or 'transcripts'"
    options = dict((k.lstrip('-'),v) for k,v in args.iteritems())

    rnacounter_main(bamname,annotname, **options)

    args['--output'].close()

