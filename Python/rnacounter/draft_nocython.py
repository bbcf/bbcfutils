#!/usr/bin/env python


"""
Usage:
   rnacounter join TAB [TAB2 ...]
   rnacounter  [...] BAM GTF
   rnacounter  [-n <int>] [-f <int>] [-s] [--nh] [--noheader] [--threshold <float>] [--gtf_ftype FTYPE]
               [--format FORMAT] [-t TYPE] [-c CHROMS] [-o OUTPUT] [-m METHOD] BAM GTF
               [--version] [-h]

Options:
   -h, --help                       Displays usage information and exits.
   -v, --version                    Displays version information and exits.
   -s, --stranded                   Compute sense and antisense reads separately [default: False].
   -n <int>, --normalize <int>      Normalization constant for RPKM. Default: (total number of mapped reads)/10^6.
   -f <int>, --fraglength <int>     Average fragment length (for transcript length correction) [default: 1].
   --nh                             Divide count by NH flag for multiply mapping reads [default: False].
   --noheader                       Remove column names from the output (helps piping) [default: False].
   --threshold <float>              Do not report counts inferior or equal to the given threshold [default: -1].
   --gtf_ftype FTYPE                Type of feature in the 3rd column of the GTF to consider [default: exon].
   --format FORMAT                  Format of the annotation file: 'gtf' or 'bed' [default: gtf].
   -t TYPE, --type TYPE             Type of genomic features to count reads in:
                                    'genes', 'transcripts', 'exons' or 'introns' [default: genes].
   -c CHROMS, --chromosomes CHROMS  Selection of chromosome names (comma-separated list).
   -o OUTPUT, --output OUTPUT       Output file to redirect stdout.
   -m METHOD, --method METHOD       Choose from 'nnls', 'raw', ('likelihood'-soon) [default: raw].
"""

import pysam
import os, sys, itertools, copy, subprocess
from operator import attrgetter
from numpy import asarray, zeros, diag, sqrt, multiply, dot
from scipy.optimize import nnls


##########################  GTF parsing  #############################


def _score(x):
    if x == '.': return 0.0
    else: return float(x)
def _strand(x):
    smap = {'+':1, '1':1, '-':-1, '-1':-1, '.':0, '0':0}
    return smap[x]

def skip_header(filename):
    """Return the number of lines starting with `#`."""
    n = 0
    with open(filename) as f:
        while f.readline()[0] == '#':
            n += 1
    return n

Ecounter = itertools.count(1)  # to give unique ids to undefined exons, see parse_gtf()
def parse_gtf(line, gtf_ftype):
    """Parse one GTF line. Return None if not an 'exon'. Return False if row is empty."""
    # GTF fields = ['chr','source','name','start','end','score','strand','frame','attributes']
    if not line: return False
    row = line.strip().split("\t")
    if len(row) < 9:
        raise ValueError("\"Attributes\" field required in GFF.")
    if row[2] != gtf_ftype:
        return None
    attrs = tuple(x.strip().split() for x in row[8].rstrip(';').split(';'))  # {gene_id: "AAA", ...}
    attrs = dict((x[0],x[1].strip("\"")) for x in attrs)
    exon_nr = Ecounter.next()
    exon_id = attrs.get('exon_id', 'E%d'%exon_nr)
    return Exon(id=(exon_nr,),
        gene_id=attrs.get('gene_id',exon_id), gene_name=attrs.get('gene_name',exon_id),
        chrom=row[0], start=max(int(row[3])-1,0), end=max(int(row[4]),0),
        name=exon_id, strand=_strand(row[6]),
        transcripts=set([attrs.get('transcript_id',exon_id)]))

def parse_bed(line, gtf_ftype):
    """Parse one BED line. Return False if *line* is empty."""
    if (not line) or (line[0]=='#') or (line[:5]=='track'): return False
    row = line.strip().split('\t')
    chrom = row[0]; start = int(row[1]); end = int(row[2]); name = row[3]
    strand = 0
    if len(row) > 4:
        if len(row) > 5: strand = _strand(row[5])
        else: strand = 0
        score = _score(row[4])
    else: score = 0.0
    exon_nr = Ecounter.next()
    #exon_id = 'E%d'%exon_nr
    return Exon(id=(exon_nr,), gene_id=name, gene_name=name, chrom=chrom, start=start, end=end,
                name=name, score=score, strand=strand, transcripts=set([name]))


##########################  Join tables  ############################


def join(tables):
    """Put 1-sample tables together into a single big table with one column per sample."""
    tabs = [open(t) for t in tables]
    out = open("joined_counts_rnacounter.txt","wb")
    lines = [t.readline().split('\t') for t in tabs]
    if any(len(L) == 0 for L in lines):
        sys.stderr.write("One of the tables is empty. Abort.")
        return 1
    try: float(lines[0][1]) # header?
    except:
        header = [lines[0][0]] + ["Count.%d"%c for c in range(len(tables))] \
                 + ["RPKM.%d"%c for c in range(len(tables))] + lines[0][3:]
        out.write('\t'.join(header))
        lines = [t.readline().split('\t') for t in tabs]
    while 1:
        if all(len(L) > 3 for L in lines):
            gid = lines[0][0]
            annot = lines[0][3:]
            cnts = [x[1] for x in lines]
            rpkm = [x[2] for x in lines]
            newline = [gid] + cnts + rpkm + annot
            out.write('\t'.join(newline))
            lines = [t.readline().split('\t') for t in tabs]
        elif all(len(L[0]) == 0 for L in lines):  # success
            break
        else:  # unequal nb of lines
            sys.stderr.write("Unequal number of lines. Abort.")
            return 1
    out.close()
    [t.close() for t in tabs]
    sys.stderr.write("Merged files into \"joined_counts_rnacounter.txt\" .\n")


#########################  Global classes  ##########################


class Counter(object):
    def __init__(self):
        self.n = 0
    def __call__(self, alignment):
        self.n += 1

class GenomicObject(object):
    def __init__(self, id=(0,),gene_id='',gene_name='',chrom='',start=0,end=0,
                 name='',score=0.0,strand=0,length=0,multiplicity=1,
                 count=0,count_anti=0,rpk=0.0,rpk_anti=0.0):
        self.id = id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand  # 1,-1,0
        self.length = length
        self.multiplicity = multiplicity
        self.count = count
        self.count_anti = count_anti
        self.rpk = rpk
        self.rpk_anti = rpk_anti
    def __and__(self,other):
        """The intersection of two GenomicObjects"""
        assert self.chrom==other.chrom, "Cannot add features from different chromosomes"
        return self.__class__(
            id = self.id + other.id,
            gene_id = '|'.join(set([self.gene_id, other.gene_id])),
            gene_name = '|'.join(set([self.gene_name, other.gene_name])),
            chrom = self.chrom,
            ##   name = '|'.join(set([self.name, other.name])),
            name = '|'.join([self.name, other.name]),
            strand = (self.strand + other.strand)/2,
            multiplicity = self.multiplicity + other.multiplicity
        )
    def __repr__(self):
        return "__%s.%s:%d-%d__" % (self.name,self.gene_name,self.start,self.end)

class Exon(GenomicObject):
    def __init__(self, transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.transcripts = transcripts   # list of transcripts it is contained in
        self.length = self.end - self.start
    def __and__(self,other):
        E = GenomicObject.__and__(self,other)
        E.transcripts = self.transcripts | other.transcripts
        return E
    def increment(self, x, alignment, multiple, stranded):
        if multiple:
            NH = 1.0
            for (k,v) in alignment.tags:
                if k=='NH':
                    NH = 1.0/v
                    break
            x = x * NH
        if stranded:
            # read/exon stand mismatch
            if (alignment.is_reverse is False and self.strand == 1) \
            or (alignment.is_reverse is True and self.strand == -1):
                self.count += x
            else:
                self.count_anti += x
        else:
            self.count += x

class Intron(Exon):
    def __init__(self, transcripts=set(), **args):
        Exon.__init__(self, **args)

class Transcript(GenomicObject):
    def __init__(self, exons=set(), **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons it contains

class Gene(GenomicObject):
    def __init__(self, exons=set(), transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons contained
        self.transcripts = transcripts   # list of transcripts contained


#####################  Operations on intervals  #####################


def intersect_exons_list(feats):
    """The intersection of a list *feats* of GenomicObjects. Pieces are unique."""
    feats = list(set(feats))
    if len(feats) == 1:
        return copy.deepcopy(feats[0])
    else:
        return reduce(lambda x,y: x&y, feats)

def cobble(exons):
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

def fuse(intervals):
    """Fuses overlapping *intervals* - a list [(a,b),(c,d),...]."""
    fused = []
    x = intervals[0]
    for y in intervals[1:]:
        if y[0] < x[1]:
            x[1] = max(x[1], y[1])
        else:
            fused.append(x)
            x = y
    fused.append(x)
    return fused

def partition_chrexons(chrexons):
    """Partition chrexons in non-overlapping chunks with distinct genes.
    The problem is that exons are sorted wrt start,end, and so the first
    exon of a gene can be separate from the second by exons of other genes
    - from the GTF we don't know how many and how far."""
    lastend = chrexons[0].end
    lastgeneids = set([chrexons[0].gene_id])
    lastindex = 0
    partition = []
    pinvgenes = {}  # map {gene_id: partitions it is found in}
    npart = 0       # partition index
    # First cut - where disjoint except if the same gene continues
    for i,exon in enumerate(chrexons):
        if (exon.start > lastend) and (exon.gene_id not in lastgeneids):
            lastend = max(exon.end,lastend)
            partition.append((lastindex,i))
            # Record in which parts the gene was found, fuse them later
            for g in lastgeneids:
                pinvgenes.setdefault(g,[]).append(npart)
            npart += 1
            lastgeneids.clear()
            lastindex = i
        else:
            lastend = max(exon.end,lastend)
        lastgeneids.add(exon.gene_id)
    partition.append((lastindex,len(chrexons)))
    for g in lastgeneids:
        pinvgenes.setdefault(g,[]).append(npart)
    # Put together intervals containing parts of the same gene mixed with others - if any
    mparts = [[p[0],p[len(p)-1]] for p in pinvgenes.itervalues() if len(p)>1]
    if mparts:
        mparts = fuse(sorted(mparts))
        toremove = set()
        for (a,b) in mparts:
            partition[b] = (partition[a][0],partition[b][1])
            toremove |= set(xrange(a,b))
        partition = [p for i,p in enumerate(partition) if i not in toremove]
    return partition

def complement(tid,tpieces):
    """From a transcript ID and its exon pieces, complement the
    intervals to get the introns of the transcript."""
    introns = []
    k = 0
    for i in range(len(tpieces)-1):
        a = tpieces[i]
        b = tpieces[i+1]
        if a.end == b.start: continue
        k += 1
        intron_id = (-1,)+a.id
        intron_name = "%s-i%d"%(tid,k)
        intron = Intron(id=intron_id, gene_id=a.gene_id, gene_name=a.gene_name, chrom=a.chrom,
            start=a.end, end=b.start, name=intron_name, strand=a.strand, transcripts=set([tid]))
        introns.append(intron)
    return introns


#############################  Counting  ############################


def toRPK(count,length,norm_cst):
    return 1000.0 * count / (length * norm_cst)
def fromRPK(rpk,length,norm_cst):
    return length * norm_cst * rpk / 1000.
def correct_fraglen_bias(rpk, length, fraglen):
    if fraglen == 1: return rpk
    newlength = max(length-fraglen+1, 1)
    return rpk * length / newlength


def count_reads(exons,ckreads,multiple,stranded):
    """Adds (#aligned nucleotides/read length) to exon counts.
    Deals with indels, junctions etc.
    :param multiple: divide the count by the NH tag.
    :param standed: for strand-specific protocols, use the strand information."""
    current_idx = 0
    nexons = len(exons)
    for alignment in ckreads:
        if current_idx >= nexons: return 0
        exon_end = exons[current_idx].end
        ali_pos = alignment.pos
        while exon_end <= ali_pos:
            current_idx += 1
            if current_idx >= nexons: return 0
            exon_end = exons[current_idx].end
        idx2 = current_idx
        exon_start = exons[idx2].start
        read_len = alignment.rlen
        ali_len = 0
        for op,shift in alignment.cigar:
            # Deletion is "from the reference", insert is "to the reference"
            if op in [0,2,3]:  # [BAM_CMATCH,BAM_CDEL,BAM_CREF_SKIP]
                # If read crosses exon left bound, trim
                if ali_pos < exon_start:
                    pos = ali_pos
                    ali_pos = min(exon_start, pos+shift)
                    shift = max(0, pos+shift-exon_start)
                # If read crosses exon right bound, maybe next exon(s)
                while ali_pos+shift >= exon_end:
                    # Score up to exon end, go to exon end, remove from shift and reset ali_len
                    if op == 0:
                        ali_len += exon_end - ali_pos
                        exons[idx2].increment(float(ali_len)/float(read_len), alignment, multiple,stranded)
                    shift -= exon_end-ali_pos
                    ali_pos = exon_end
                    ali_len = 0
                    # Next exon
                    idx2 += 1
                    if idx2 >= nexons: return 0
                    exon_start = exons[idx2].start
                    exon_end = exons[idx2].end
                    # If op crosses exon left bound, go to exon start with rest of shift
                    # If op ends before the next exon, go to end of op and reset shift
                    if ali_pos < exon_start:
                        pos = ali_pos
                        ali_pos = min(exon_start, pos+shift)
                        shift = max(0, pos+shift-exon_start)
                # If a bit of op remains overlapping the next exon
                if op == 0:
                    ali_len += shift
                ali_pos += shift  # got to start of next op in prevision for next round
            elif op == 1:  # BAM_CINS
                ali_len += shift;
            exons[idx2].increment(float(ali_len)/float(read_len), alignment, multiple,stranded)


def get_total_nreads(sam):
    Ncounter = Counter()
    for ref,length in itertools.izip(sam.references,sam.lengths):
        sam.fetch(ref,0,length, callback=Ncounter)
    return Ncounter.n


######################  Expression inference  #######################


def is_in(feat_class,x,feat_id):
    """Returns True if Exon *x* is part of the gene/transcript *feat_id*."""
    if feat_class == Transcript:
        return feat_id in x.transcripts
    elif feat_class == Gene:
        return feat_id in x.gene_id.split('|')
    else:
        return x.name in feat_id.split('|') or x.name == feat_id
            # x is an exon: x == itself or x contains the piece
            # x is a piece: p == feat_id

def estimate_expression_NNLS(feat_class, pieces, ids, exons, norm_cst, stranded, weighted):
    #--- Build the exons-transcripts structure matrix:
    # Lines are exons, columns are transcripts,
    # so that A[i,j]!=0 means "transcript Tj contains exon Ei".
    n = len(pieces)
    m = len(ids)
    A = zeros((n,m))
    for i,p in enumerate(pieces):
        for j,f in enumerate(ids):
            A[i,j] = 1. if is_in(feat_class,p,f) else 0.
    #--- Build the exons RPKs vector
    E = asarray([p.rpk for p in pieces])
    if weighted:
        w = sqrt(asarray([p.length for p in pieces]))
        W = diag(w)
        A = dot(W,A)
        E = multiply(E, w)
    #--- Solve for transcripts RPK
    T,rnorm = nnls(A,E)
    #-- Same for antisense if stranded protocol
    if stranded:
        E_anti = asarray([p.rpk_anti for p in pieces])
        T_anti,rnorm_anti = nnls(A,E_anti)
    #--- Store result in *feat_class* objects
    feats = []
    frpk_anti = fcount_anti = 0.0
    for i,f in enumerate(ids):
        exs = sorted([e for e in exons if is_in(feat_class,e,f)], key=attrgetter('start','end'))
        flen = sum([p.length for p in pieces if is_in(feat_class,p,f)])
        frpk = T[i]
        fcount = fromRPK(T[i],flen,norm_cst)
        if stranded:
            frpk_anti = T_anti[i]
            fcount_anti = fromRPK(T_anti[i],flen,norm_cst)
        feats.append(feat_class(name=f, length=flen,
                rpk=frpk, rpk_anti=frpk_anti, count=fcount, count_anti=fcount_anti,
                chrom=exs[0].chrom, start=exs[0].start, end=exs[len(exs)-1].end,
                gene_id=exs[0].gene_id, gene_name=exs[0].gene_name, strand=exs[0].strand))
    return feats


def estimate_expression_raw(feat_class, pieces, ids, exons, norm_cst, stranded):
    """For each feature ID in *ids*, just sum the score of its components as one
    commonly does for genes from exon counts."""
    feats = []
    frpk_anti = fcount_anti = 0.0
    for i,f in enumerate(ids):
        exs = sorted([e for e in exons if is_in(feat_class,e,f)], key=attrgetter('start','end'))
        inner = [p for p in pieces if (len(p.gene_id.split('|'))==1 and is_in(feat_class,p,f))]
        if len(inner)==0:
            flen = 1
            fcount = frpk = 0.0
        else:
            flen = sum([p.length for p in inner])
            fcount = sum([p.count for p in inner])
            frpk = toRPK(fcount,flen,norm_cst)
            if stranded:
                fcount_anti = sum([p.count_anti for p in inner])
                frpk_anti = toRPK(fcount_anti,flen,norm_cst)
        feats.append(feat_class(name=f, length=flen,
                rpk=frpk, rpk_anti=frpk_anti, count=fcount, count_anti=fcount_anti,
                chrom=exs[0].chrom, start=exs[0].start, end=exs[len(exs)-1].end,
                gene_id=exs[0].gene_id, gene_name=exs[0].gene_name, strand=exs[0].strand))
    return feats


###########################  Main script  ###########################


def filter_transcripts(t2p, readlength):
    """*t2p* is a map {transcriptID: [exon pieces]}.
    Find transcripts that differ from others by less than a few exons of less than
    a read length. Return the set of IDs of redundant transcripts."""
    seen = set()  # transcript structures, as tuples of exon ids
    toremove = set() # too close transcripts
    for t,tpieces in t2p.iteritems():
        filtered_ids = tuple([tp.id for tp in tpieces if tp.length < readlength])
        if filtered_ids in seen:
            toremove.add(t)
        else:
            seen.add(filtered_ids)
    return toremove


def process_chunk(ckexons, sam, chrom, options):
    """Distribute counts across transcripts and genes of a chunk *ckexons*
    of non-overlapping exons."""

    norm_cst = options['normalize']
    stranded = options['stranded']
    output = options['output']
    types = options['type']
    methods = options['method']
    threshold = options['threshold']
    fraglength = options['fraglength']
    readlength = options["readlength"]
    weighted = True

    #--- Regroup occurrences of the same Exon from a different transcript
    exons = []
    for key,group in itertools.groupby(ckexons, attrgetter('id')):
        # ckexons are sorted by id because chrexons were sorted by chrom,start,end
        exon0 = group.next()
        for gr in group:
            exon0.transcripts.add(gr.transcripts.pop())
        exons.append(exon0)
    gene_ids = set(e.gene_id for e in exons)
    exon_names = set(e.name for e in exons)

    #--- Cobble all these intervals
    pieces = cobble(exons)  # sorted

    #--- Build transcript to pieces mapping
    t2p = {}
    transcript_ids = set()
    for p in pieces:
        for t in p.transcripts:
            transcript_ids.add(t)
            t2p.setdefault(t,[]).append(p)

    #--- Filter out too similar transcripts
    if 1 in types or 3 in types:
        toremove = filter_transcripts(t2p, readlength)
        for t in toremove:
            transcript_ids.remove(t)
            t2p.pop(t)
        for p in pieces:
            p.transcripts = set(t for t in p.transcripts if t not in toremove)

    #--- Count reads in each piece
    lastend = max(e.end for e in exons)
    ckreads = sam.fetch(chrom, exons[0].start, lastend)
    count_reads(pieces,ckreads,options['nh'],stranded)

    #--- Same for introns, if selected
    intron_pieces = []
    if 3 in types:
        introns = []
        for tid,tpieces in sorted(t2p.items()):
            introns.extend(complement(tid,tpieces))  # tpieces is already sorted
        if introns:
            intron_exon_pieces = cobble(introns+exons)
            intron_pieces = [ip for ip in intron_exon_pieces \
                             if ip.length > readlength \
                             and not any([n in exon_names for n in ip.name.split('|')])]
            if intron_pieces:
                lastend = max([intron.end for intron in introns])
                ckreads = sam.fetch(chrom, intron_pieces[0].start, lastend)
                count_reads(intron_pieces,ckreads,options['nh'],stranded)
        for ip in intron_pieces:   # transcripts is a set, so intron names join randomly...
            ip.name = '|'.join(sorted(ip.name.split('|')))

    #--- Calculate RPK
    for p in itertools.chain(pieces,intron_pieces):
        p.rpk = toRPK(p.count,p.length,norm_cst)
    if stranded:
        for p in itertools.chain(pieces,intron_pieces):
            p.rpk_anti = toRPK(p.count_anti,p.length,norm_cst)

    #--- Infer gene/transcript counts
    for p in pieces:
        p.name = '|'.join(sorted(set(p.name.split('|'))))  # remove duplicates in names
    genes=[]; transcripts=[]; exons2=[]; introns2=[]
    # Genes - 0
    if 0 in types:
        method = methods.get(0,0)
        if method == 0:
            genes = estimate_expression_raw(Gene,pieces,gene_ids,exons,norm_cst,stranded)
        elif method == 1:
            genes = estimate_expression_NNLS(Gene,pieces,gene_ids,exons,norm_cst,stranded,weighted)
        for gene in genes:
            gene.rpk = correct_fraglen_bias(gene.rpk, gene.length, fraglength)
    # Transcripts - 1
    if 1 in types:
        method = methods.get(1,0)
        if method == 1:
            transcripts = estimate_expression_NNLS(Transcript,pieces,transcript_ids,exons,norm_cst,stranded,weighted)
        elif method == 0:
            transcripts = estimate_expression_raw(Transcript,pieces,transcript_ids,exons,norm_cst,stranded)
        for trans in transcripts:
            trans.rpk = correct_fraglen_bias(trans.rpk, trans.length, fraglength)
    # Exons - 2
    if 2 in types:
        method = methods.get(2,0)
        if method == 0:
            exons2 = list(pieces)   # !
        elif method == 1:
            exon_ids = set([e.name for e in exons])
            exons2 = estimate_expression_NNLS(Exon,pieces,exon_ids,exons,norm_cst,stranded)
    # Introns - 3
    if 3 in types and intron_pieces:
        method = methods.get(3,0)
        if method == 0:
            introns2 = list(intron_pieces)   # !
        elif method == 1:
            intron_ids = set(['|'.join(sorted(e.name.split('|'))) for e in introns])
            introns2 = estimate_expression_NNLS(Exon,intron_pieces,intron_ids,introns,norm_cst,stranded)

    print_output(output, genes,transcripts,exons2,introns2, threshold,stranded)


def print_output(output, genes,transcripts,exons,introns, threshold,stranded):
    igenes = itertools.ifilter(lambda x:x.count > threshold, genes)
    itranscripts = itertools.ifilter(lambda x:x.count > threshold, transcripts)
    iexons = itertools.ifilter(lambda x:x.count > threshold, exons)
    iintrons = itertools.ifilter(lambda x:x.count > threshold, introns)
    if stranded:
        for f in itertools.chain(igenes,itranscripts,iexons,iintrons):
            towrite = [str(x) for x in [f.name,f.count,f.rpk,f.chrom,f.start,f.end,
                                        f.strand,f.gene_name,f.__class__.__name__.lower(),'sense']]
            output.write('\t'.join(towrite)+'\n')
            towrite = [str(x) for x in [f.name,f.count_anti,f.rpk_anti,f.chrom,f.start,f.end,
                                        f.strand,f.gene_name,f.__class__.__name__.lower(),'antisense']]
            output.write('\t'.join(towrite)+'\n')
    else:
        for f in itertools.chain(igenes,itranscripts,iexons,iintrons):
            towrite = [str(x) for x in [f.name,f.count,f.rpk,f.chrom,f.start,f.end,
                                        f.strand,f.gene_name,f.__class__.__name__.lower()]]
            output.write('\t'.join(towrite)+'\n')


def rnacounter_main(bamname, annotname, options):
    # Index BAM if necessary
    if not os.path.exists(bamname+'.bai'):
        sys.stderr.write("BAM index not found. Indexing...")
        subprocess.check_call("samtools index %s" % bamname, shell=True)
        sys.stderr.write("...done.\n")

    # Open SAM. Get read length
    sam = pysam.Samfile(bamname, "rb")
    options["readlength"] = sam.next().rlen
    sam.close()
    sam = pysam.Samfile(bamname, "rb")

    # Open GTF. Skip header lines
    nhead = skip_header(annotname)
    annot = open(annotname, "r")
    for _ in range(nhead): annot.readline()

    if options['output'] is None: options['output'] = sys.stdout
    else: options['output'] = open(options['output'], "wb")
    if options['noheader'] is False:
        header = ['ID','Count','RPKM','Chrom','Start','End','Strand','GeneName','Type']
        if options['stranded']: header += ['Sense']
        options['output'].write('\t'.join(header)+'\n')

    if options['format'] == 'gtf':
        parse = parse_gtf
    elif options['format'] == 'bed':
        parse = parse_bed

    # Cross 'chromosomes' option with available BAM headers
    if options['chromosomes']:
        chromosomes = [c for c in sam.references if c in options['chromosomes']]
    else:
        chromosomes = sam.references

    ## Get total number of reads
    if options['normalize'] is None:
        options['normalize'] = get_total_nreads(sam) / 1.0e6
    else:
        options['normalize'] = float(options['normalize'])
    options['normalize'] = 1.0

    # Initialize
    gtf_ftype = options['gtf_ftype']
    chrom = ''
    while chrom not in chromosomes:
        exon = None
        while exon is None:
            row = annot.readline().strip()
            exon = parse(row, gtf_ftype)
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
                exon = parse(row, gtf_ftype)  # None if not an exon, False if EOF
            if not row: break
            chrom = exon.chrom
        if chrexons:
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
    assert os.path.exists(bamname), "BAM file not found: %s" %bamname
    assert os.path.exists(annotname), "GTF file not found: %s" %annotname

    if args['--chromosomes'] is None: args['--chromosomes'] = []
    else: args['--chromosomes'] = args['--chromosomes'].split(',')

    assert args['--format'].lower() in ['gtf','bed'], \
        "FORMAT must be one of 'gtf' or 'bed'."

    # Type: one can actually give both as "-t genes,transcripts" but they
    # will be mixed in the output stream. Split the output using the last field ("Type").
    args['--type'] = [x.lower() for x in args['--type'].split(',')]
    assert all(x in ["genes","transcripts","exons","introns"] for x in args['--type']), \
        "TYPE must be one of 'genes', 'transcripts', 'exons' or 'introns'."
    type_map = {'genes':0, 'transcripts':1, 'exons':2, 'introns':3}  # avoid comparing strings later
    args['--type'] = [type_map[x] for x in args['--type']]

    # Same for methods. If given as a list, the length must be that of `--type`,
    # the method at index i will be applied to feature type at index i.
    args['--method'] = [x.lower() for x in args['--method'].split(',')]
    if len(args['--method']) > 1:  # multiple methods given
        assert len(args['--method']) == len(args['--type']), \
            "TYPE and METHOD arguments must have the same number of elements."
    elif len(args['--type']) > 1:  # apply same method to all types
        args['--method'] = args['--method'] * len(args['--type'])
    assert all(x in ["raw","nnls","wnnls"] for x in args['--method']), \
        "METHOD must be one of 'raw' or 'nnls'."
    method_map = {'raw':0, 'nnls':1, "wnnls":2}  # avoid comparing strings later
    args['--method'] = [method_map[x] for x in args['--method']]
    args['--method'] = dict(zip(args['--type'],args['--method']))

    try: args['--threshold'] = float(args['--threshold'])
    except ValueError: raise ValueError("--threshold must be numeric.")
    try: args['--fraglength'] = int(args['--fraglength'])
    except ValueError: raise ValueError("--fraglength must be an integer.")

    options = dict((k.lstrip('-').lower(), v) for k,v in args.iteritems())
    return bamname, annotname, options


if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    if args['join']:
        join([args['TAB']]+args['TAB2'])
    else:
        bamname, annotname, options = parse_args(args)
        rnacounter_main(bamname,annotname, options)


#----------------------------------------------#
# This code was written by Julien Delafontaine #
# EPFL,BBCF: http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch                       #
#----------------------------------------------#
