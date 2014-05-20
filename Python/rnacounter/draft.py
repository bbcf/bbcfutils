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

import pysam
import os, sys
import itertools
from numpy import asarray, zeros
import copy
from scipy.optimize import nnls


_TEST_ = False


Ecounter = itertools.count(1)  # to give unique ids to undefined exons, see parse_gtf()
def parse_gtf(row):
    """Parse one GTF line. Return None if not an 'exon'. Return False if row is empty."""
    # GTF fields = ['chr','source','name','start','end','score','strand','frame','attributes']
    def _score(x):
        if str(x) == '.': return 0
        else: return float(x)
    def _strand(x):
        smap = {'+':1, 1:1, '-':-1, -1:-1, '.':0, 0:0}
        return smap[x]
    if not row: return False
    row = row.strip().split("\t")
    if len(row) < 9:
        raise ValueError("\"Attributes\" field required in GFF.")
    if row[2] != 'exon':
        return None
    attrs = (x.strip().split() for x in row[8].split(';'))  # {gene_id: "AAA", ...}
    attrs = dict((x[0],x[1].strip("\"")) for x in attrs)
    exon_id = attrs.get('exon_id', 'E%d'%Ecounter.next())
    return Exon(id=exon_id, gene_id=attrs['gene_id'], gene_name=attrs['gene_name'],
                chrom=row[0], start=int(row[3])-1, end=int(row[4]),
                name=exon_id, score=_score(row[5]), strand=_strand(row[6]),
                transcripts=[attrs['transcript_id']])


class GenomicObject(object):
    def __init__(self, id='',gene_id='',gene_name='',chrom='',start=0,end=0,
                 name='',score=0.0,count=0,count_rev=0,rpk=0.0,strand=0,length=0,seq='',multiplicity=1):
        self.id = id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        #self.score = score
        self.count = count
        self.count_rev = count_rev
        self.rpk = rpk
        self.strand = strand  # 1,-1,0
        self.length = length
        #self.seq = seq  # sequence
        self.multiplicity = multiplicity
    def __and__(self,other):
        """The intersection of two GenomicObjects"""
        assert self.chrom==other.chrom, "Cannot add features from different chromosomes"
        selfid = (self.id,) if isinstance(self.id,int) else self.id
        otherid = (other.id,) if isinstance(other.id,int) else other.id
        return self.__class__(
            id = selfid + otherid,
            gene_id = '|'.join(set([self.gene_id, other.gene_id])),
            gene_name = '|'.join(set([self.gene_name, other.gene_name])),
            chrom = self.chrom,
            #start = max(self.start, other.start),
            #end = min(self.end, other.end),
            ##   name = '|'.join(set([self.name, other.name])),
            name = '|'.join([self.name, other.name]),
            #score = self.score + other.score,
            strand = (self.strand + other.strand)/2,
            #length = min(self.end, other.end) - max(self.start, other.start),
            multiplicity = self.multiplicity + other.multiplicity
        )
    def __repr__(self):
        return "<%s (%d-%d) %s>" % (self.name,self.start,self.end,self.gene_name)

class Exon(GenomicObject):
    def __init__(self, transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.transcripts = transcripts   # list of transcripts it is contained in
        self.length = self.end - self.start
    def __and__(self,other):
        E = GenomicObject.__and__(self,other)
        E.transcripts = set(self.transcripts) | set(other.transcripts)
        return E
    def increment(self, x, alignment, multiple=False, stranded=False):
        if multiple:
            NH = [1.0/t[1] for t in alignment.tags if t[0]=='NH']+[1]
            x = x*NH[0]
        if stranded:
            # read/exon stand mismatch
            if alignment.is_reverse and self.strand != -1:
                self.count_rev += x
            else:
                self.count += x
        else:
            self.count += x


class Transcript(GenomicObject):
    def __init__(self, exons=[], **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons it contains

class Gene(GenomicObject):
    def __init__(self, exons=set(),transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons contained
        self.transcripts = transcripts   # list of transcripts contained


def intersect_exons_list(feats, multiple=False):
    """The intersection of a list *feats* of GenomicObjects.
    If *multiple* is True, permits multiplicity: if the same exon E1 is
    given twice, there will be "E1|E1" parts. Otherwise pieces are unique."""
    if multiple is False:
        feats = list(set(feats))
    if len(feats) == 1:
        return copy.deepcopy(feats[0])
    else:
        return reduce(lambda x,y: x&y, feats)


def cobble(exons, multiple=False):
    """Split exons into non-overlapping parts.
    :param multiple: see intersect_exons_list()."""
    ends = [(e.start,1,e) for e in exons] + [(e.end,0,e) for e in exons]
    ends.sort()
    active_exons = []
    cobbled = []
    for i in xrange(len(ends)-1):
        a = ends[i]
        b = ends[i+1]
        if a[1]==1:
            active_exons.append(a[2])
        elif a[1]==0:
            active_exons.remove(a[2])
        if len(active_exons)==0:
            continue
        if a[0]==b[0]:
            continue
        e = intersect_exons_list(active_exons)
        e.start = a[0]; e.end = b[0]; e.length = b[0]-a[0]
        cobbled.append(e)
    return cobbled


######################################################################


def process_chrexons(chrexons, sam,chrom, multiple, stranded):
    # Process chunks of overlapping exons / exons of the same gene
    lastend = chrexons[0].end
    lastgeneid = ''
    ckexons = []
    for exon in chrexons:
        # Store
        if (exon.start <= lastend) or (exon.gene_id == lastgeneid):
            ckexons.append(exon)
        # Process the stored chunk of exons
        else:
            ckgenes,cktranscripts = process_chunk(ckexons, sam, chrom, lastend, multiple, stranded)
            ckexons = [exon]
        lastend = max(exon.end,lastend)
        lastgeneid = exon.gene_id
    ckgenes,cktranscripts = process_chunk(ckexons, sam, chrom, lastend, multiple,stranded)


def process_chunk(ckexons, sam, chrom, lastend, multiple, stranded):
    """Distribute counts across transcripts and genes of a chunk *ckexons*
    of non-overlapping exons."""

    def toRPK(count,length):
        return 1000.0 * count / length
    def fromRPK(rpk,length):
        return length * rpk / 1000.


    #--- Regroup occurrences of the same Exon from a different transcript
    exons = []
    for key,group in itertools.groupby(ckexons, lambda x:x.id):
        # ckexons are sorted because chrexons were sorted by chrom,start,end
        exon0 = group.next()
        for g in group:
            exon0.transcripts.append(g.transcripts[0])
        exons.append(exon0)
    del ckexons


    #--- Cobble all these intervals
    pieces = cobble(exons)


    #--- Filter out too similar transcripts,
    # e.g. made of the same exons up to 100bp.
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


    #--- Get all reads from this chunk - iterator
    ckreads = sam.fetch(chrom, exons[0].start, lastend)


    #--- Count reads in each piece -- from rnacounter.cc
    def count_reads(exons,ckreads,multiple=False,stranded=False):
        """Adds (#aligned nucleotides/read length) to exon counts.
        Deals with indels, junctions etc.
        :param multiple: divide the count by the NH tag.
        :param standed: for strand-specific protocols, use the strand information."""
        current_pos = 0   # exon_list.current_pos -- must be a pointer
        for alignment in ckreads:
            if current_pos > len(exons): return 0
            exon_end = exons[current_pos].end
            ali_pos = alignment.pos
            while exon_end <= ali_pos:
                current_pos += 1
                if current_pos > len(exons): return 0
                exon_end = exons[current_pos].end
            pos2 = current_pos
            exon_start = exons[pos2].start
            read_len = alignment.rlen
            ali_len = 0
            for op,shift in alignment.cigar:
                if op in [0,2,3]:  # [BAM_CMATCH,BAM_CDEL,BAM_CREF_SKIP]
                    # If read crosses exon left bound
                    if ali_pos < exon_start:
                        ali_pos = min(exon_start, ali_pos+shift)
                        shift = max(0, ali_pos+shift-exon_start)  # part of the read overlapping
                    # If read crosses exon right bound, maybe next exon(s)
                    while ali_pos+shift >= exon_end:
                        if op == 0: ali_len += exon_end - ali_pos
                        exons[pos2].increment(float(ali_len)/float(read_len), alignment, multiple, stranded)
                        shift -= exon_end-ali_pos  # remaining part of the read
                        ali_pos = exon_end
                        ali_len = 0
                        pos2 += 1
                        if pos2 >= len(exons): return 0
                        exon_start = exons[pos2].start
                        exon_end = exons[pos2].end
                        if ali_pos < exon_start:
                            ali_pos = min(exon_start, ali_pos+shift)
                            shift = max(0, ali_pos+shift-exon_start)  # from exon start to end of the read
                    ali_pos += shift
                    if op == 0: ali_len += shift
                elif op == 1:  # BAM_CINS
                    ali_len += shift;
            if ali_len < 1: return 0;
            exons[pos2].increment(float(ali_len)/float(read_len), alignment, multiple, stranded)

    count_reads(pieces,ckreads,multiple,stranded)
    #--- Calculate RPK
    for p in pieces:
        p.rpk = toRPK(p.count,p.length)


    def estimate_expression(feat_class, pieces, ids):
        #--- Build the exons-transcripts structure matrix:
        # Lines are exons, columns are transcripts,
        # so that A[i,j]!=0 means "transcript Tj contains exon Ei".
        if feat_class == Gene:
            is_in = lambda p,g: g in p.gene_id.split('|')
        elif feat_class == Transcript:
            is_in = lambda p,t: t in p.transcripts
        n = len(pieces)
        m = len(ids)
        A = zeros((n,m))
        for i,p in enumerate(pieces):
            for j,f in enumerate(ids):
                A[i,j] = 1 if is_in(p,f) else 0
        #--- Build the exons scores vector
        E = asarray([p.rpk for p in pieces])
        #--- Solve for RPK
        T,rnorm = nnls(A,E)
        #--- Store result in *feat_class* objects
        feats = []
        for i,f in enumerate(ids):
            exs = sorted([p for p in pieces if is_in(p,f)], key=lambda x:(x.start,x.end))
            flen = sum(p.length for p in pieces if is_in(p,f))
            feats.append(Transcript(name=f, start=exs[0].start, end=exs[-1].end,
                    length=flen, rpk=T[i], count=fromRPK(T[i],flen),
                    chrom=exs[0].chrom, gene_id=exs[0].gene_id, gene_name=exs[0].gene_name))
        return feats

    genes = estimate_expression(Gene, pieces, gene_ids)
    transcripts = estimate_expression(Transcript, pieces, transcript_ids)

    # Test
    if _TEST_ and exons[0].gene_name == "Gapdh":
        print "Transcripts:"
        for t in transcripts: print t,t.count,t.rpk
        print "Gene:"
        for g in genes: print g,g.count,g.rpk
        print "Pieces:"
        for p in pieces:print p,p.rpk

    return genes,transcripts


def rnacounter_main(bamname, annotname, multiple=False, stranded=False, output=sys.stdout,
                  normalize=False, chromosomes=[], fraglength=0):

    sam = pysam.Samfile(bamname, "rb")
    annot = open(annotname, "r")

    if len(chromosomes) > 0: chromosomes = [c for c in sam.references if c in chromosomes]
    else: chromosomes = sam.references

    chrom = ''
    while chrom not in chromosomes:
        exon = None
        while exon is None:
            row = annot.readline().strip()
            exon = parse_gtf(row)
        chrom = exon.chrom
    lastchrom = chrom

    while row:
        # Load all GTF exons of a chromosome in memory, sort and process
        chrexons = []
        print ">> Chromosome", chrom
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
            process_chrexons(chrexons,sam,lastchrom, multiple,stranded)
        lastchrom = chrom

    annot.close()
    sam.close


######################################################################


from docopt import docopt

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




# Gapdh id: ENSMUSG00000057666
# Gapdh transcripts: ENSMUST00000147954, ENSMUST00000147954, ENSMUST00000118875
#                    ENSMUST00000073605, ENSMUST00000144205, ENSMUST00000144588

