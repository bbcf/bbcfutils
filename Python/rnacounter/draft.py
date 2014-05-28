"""
Count reads on genes and transcripts from a genome-level BAM file and a
GTF/GFF file describing the exons, such as those provided by Emsembl or GenRep.
The GTF is assumed to be sorted at least w.r.t. chromosome name,
and the chromosome identifiers in the GTF must be the same as the BAM references.

If your GTF does not represent exons but custom genomic intervals to simply count
reads in, provide at least a unique `exon_id` in the attributes as a feature name,
and the type field (column 2) must be set to 'exon'.
If not specified, `gene_id`, `transcript_id` and `exon_id` will all get the value
of `exon_id` and be considered as independant features.

Usage:
   rnacounter  [-t TYPE] [-n <int>] [-l <int>] [-s] [-m] [-c CHROMS] [-o OUTPUT] BAM GTF
               [--version] [-h]

Options:
   -t TYPE, --type TYPE             Type of genomic feature to count on: 'genes' or 'transcripts' [default: genes].
   -n <int>, --normalize <int>      Normalization constant for RPKM. Default: (total number of mapped reads)/10^6.
   -l <int>, --fraglength <int>     Average fragment length [default: 350].
   -s, --stranded                   Compute sense and antisense reads separately [default: False].
   -m, --multiple                   Divide count by NH flag for multiply mapping reads [default: False].
   -c CHROMS, --chromosomes CHROMS  Selection of chromosome names (comma-separated list).
   -o OUTPUT, --output OUTPUT       Output file to redirect stdout.
   -v, --version                    Displays version information and exits.
   -h, --help                       Displays usage information and exits.
"""

import pysam
import os, sys, itertools

from rnacounter0 import count_reads, partition_chrexons, estimate_expression, get_total_nreads
from rnacounter0 import _score, _strand, toRPK, cobble
from rnacounter0 import Exon, Gene, Transcript


Ecounter = itertools.count(1)  # to give unique ids to undefined exons, see parse_gtf()
def parse_gtf(row):
    """Parse one GTF line. Return None if not an 'exon'. Return False if row is empty."""
    # GTF fields = ['chr','source','name','start','end','score','strand','frame','attributes']
    if not row: return False
    row = row.strip().split("\t")
    if len(row) < 9:
        raise ValueError("\"Attributes\" field required in GFF.")
    if row[2] != 'exon':
        return None
    attrs = (x.strip().split() for x in row[8].split(';'))  # {gene_id: "AAA", ...}
    attrs = dict((x[0],x[1].strip("\"")) for x in attrs)
    exon_id = attrs.get('exon_id', 'E%d'%Ecounter.next())
    return Exon(id=exon_id,
        gene_id=attrs.get('gene_id',exon_id), gene_name=attrs.get('gene_name',exon_id),
        chrom=row[0], start=int(row[3])-1, end=int(row[4]),
        name=exon_id, score=_score(row[5]), strand=_strand(row[6]),
        transcripts=[attrs.get('transcript_id',exon_id)], exon_number=int(attrs.get('exon_number',1)))


######################################################################


def process_chunk(ckexons, sam, chrom, options):
    """Distribute counts across transcripts and genes of a chunk *ckexons*
    of non-overlapping exons."""

    #--- Regroup occurrences of the same Exon from a different transcript
    exons = []
    for key,group in itertools.groupby(ckexons, lambda x:x.id):
        # ckexons are sorted by id because chrexons were sorted by chrom,start,end
        exon0 = group.next()
        for g in group:
            exon0.transcripts.append(g.transcripts[0])
        exons.append(exon0)
    del ckexons

    #--- Cobble all these intervals
    pieces = cobble(exons)

    #--- Filter out too similar transcripts, e.g. made of the same exons up to 100bp.
    t2e = {}                               # map {transcript: [pieces IDs]}
    for p in pieces:
        if p.length < 100: continue        # filter out cobbled pieces of less that read length
        for t in p.transcripts:
            t2e.setdefault(t,[]).append(p.id)
    e2t = {}
    for t,e in t2e.iteritems():
        es = tuple(sorted(e))              # combination of pieces indices
        e2t.setdefault(es,[]).append(t)    # {(pieces IDs combination): [transcripts with same struct]}
    # Replace too similar transcripts by the first of the list, arbitrarily
    transcript_ids = set()  # full list of remaining transcripts
    tx_replace = dict((badt,tlist[0]) for tlist in e2t.values() for badt in tlist[1:] if len(tlist)>1)
    for p in pieces:
        filtered = set([tx_replace.get(t,t) for t in p.transcripts])
        transcript_ids |= filtered
        p.transcripts = list(filtered)
    transcript_ids = list(transcript_ids)
    gene_ids = list(set(e.gene_id for e in exons))

    #print ">> Process chunk with genes:", gene_ids

    #--- Get all reads from this chunk - iterator
    lastend = max(e.end for e in exons)
    ckreads = sam.fetch(chrom, exons[0].start, lastend)

    #--- Count reads in each piece -- from rnacounter.cc
    count_reads(pieces,ckreads,options['multiple'],options['stranded'])

    #--- Calculate RPK
    for p in pieces:
        p.rpk = toRPK(p.count,p.length,options['normalize'])

    #--- Print output
    genes = []; transcripts = []
    if 'genes' in options['type']:
        genes = estimate_expression(Gene, pieces, gene_ids, exons, options['normalize'])
    if 'transcripts' in options['type']:
        transcripts = estimate_expression(Transcript, pieces, transcript_ids, exons, options['normalize'])
    for f in itertools.chain(genes,transcripts):
        towrite = [str(x) for x in [f.name,f.count,f.rpk,f.chrom,f.start,f.end,
                                    f.strand,f.gene_name,f.__class__.__name__.lower()]]
        options['output'].write('\t'.join(towrite)+'\n')


def rnacounter_main(bamname, annotname, options):
    sam = pysam.Samfile(bamname, "rb")
    annot = open(annotname, "r")

    if options['output'] is None: options['output'] = sys.stdout
    else: options['output'] = open(options['output'], "wb")

    # Cross 'chromosomes' option with available BAM headers
    if len(options['chromosomes']) > 0:
        chromosomes = [c for c in sam.references if c in options['chromosomes']]
    else:
        chromosomes = sam.references

    # Get total number of reads
    if options['normalize'] is None:
        options['normalize'] = get_total_nreads(sam) / 1.0e6
    else:
        options['normalize'] = float(options['normalize'])

    # Initialize
    chrom = ''
    while chrom not in chromosomes:
        exon = None
        while exon is None:
            row = annot.readline().strip()
            exon = parse_gtf(row)
        chrom = exon.chrom
    lastchrom = chrom

    # Process together all exons of one chromosome at a time
    while row:
        chrexons = []
        while chrom == lastchrom:
            if (exon.end - exon.start > 1) and (exon.chrom in chromosomes):
                chrexons.append(exon)
            # Fetch next exon
            exon = None
            while exon is None:
                row = annot.readline().strip()
                exon = parse_gtf(row)  # None if not an exon, False if EOF
            if not row: break
            chrom = exon.chrom
        if len(chrexons) > 0:
            chrexons.sort(key=lambda x: (x.start,x.end))
            partition = partition_chrexons(chrexons)
            # Process chunks
            for (a,b) in partition:
                process_chunk(chrexons[a:b], sam, lastchrom, options)
        lastchrom = chrom

    options['output'].close()
    annot.close()
    sam.close()


######################################################################


from docopt import docopt

def usage_string():
    return __doc__

def parse_args(args):
    bamname = os.path.abspath(args['BAM'])
    annotname = os.path.abspath(args['GTF'])
    assert(os.path.exists(bamname))
    assert(os.path.exists(annotname))

    if args['--chromosomes'] is None: args['--chromosomes'] = []
    else: args['--chromosomes'] = args['--chromosomes'].split(',')

    # Type: one can actually give both as "-t genes,transcripts" but they
    # will be mixed in the output stream. Split the output using the last field ("Type").
    args['--type'] = [x.lower() for x in args['--type'].split(',')]
    assert all(x in ["genes","transcripts"] for x in args['--type']), \
        "TYPE must be one of 'genes' or 'transcripts'"

    options = dict((k.lstrip('-'),v) for k,v in args.iteritems())
    return bamname, annotname, options


if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    bamname, annotname, options = parse_args(args)
    rnacounter_main(bamname,annotname, options)


#----------------------------------------------#
# This code was written by Julien Delafontaine #
# EPFL,BBCF: http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch                       #
#----------------------------------------------#
