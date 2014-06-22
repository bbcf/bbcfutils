#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True
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
   rnacounter  [-t TYPE] [-n <int>] [-l <int>] [-s] [-m] [-c CHROMS] [-o OUTPUT]
               [--format FORMAT] [--method METHOD] BAM GTF
               [--version] [-h]

Options:
   -f FORMAT, --format FORMAT       Format of the annotation file: 'gtf' or 'bed' [default: gtf].
   -t TYPE, --type TYPE             Type of genomic feature to count on: 'genes' or 'transcripts' [default: genes].
   -n <int>, --normalize <int>      Normalization constant for RPKM. Default: (total number of mapped reads)/10^6.
   -l <int>, --fraglength <int>     Average fragment length [default: 350].
   -s, --stranded                   Compute sense and antisense reads separately [default: False].
   -m, --multiple                   Divide count by NH flag for multiply mapping reads [default: False].
   -c CHROMS, --chromosomes CHROMS  Selection of chromosome names (comma-separated list).
   -o OUTPUT, --output OUTPUT       Output file to redirect stdout.
   --method METHOD                  Choose from 'nnls', 'raw', ('likelihood'-soon) [default: raw].
   -v, --version                    Displays version information and exits.
   -h, --help                       Displays usage information and exits.
"""

import pysam
import os, sys, itertools, copy
from operator import attrgetter
from numpy import asarray, zeros
from scipy.optimize import nnls

import numpy as np
cimport numpy as cnp
DTYPE = np.double               # fix a datatype for the arrays
ctypedef cnp.double_t DTYPE_t   # assign a corresponding compile-time C type to DTYPE_t

cdef extern from "math.h":
    double sqrt(double x)
    double pow(double x,double y)


##########################  GTF parsing  #############################


cdef inline double _score(x):
    if x == '.': return 0.0
    else: return float(x)
cdef inline int _strand(x):
    smap = {'+':1, 1:1, '-':-1, -1:-1, '.':0, 0:0}
    return smap[x]
Ecounter = itertools.count(1)  # to give unique ids to undefined exons, see parse_gtf()

def parse_gtf(str line):
    """Parse one GTF line. Return None if not an 'exon'. Return False if *line* is empty."""
    # GTF fields = ['chr','source','name','start','end','score','strand','frame','attributes']
    cdef list row
    if (not line): return False
    row = line.strip().split("\t")
    if len(row) < 9:
        raise ValueError("\"Attributes\" field required in GFF.")
    if row[2] != 'exon':
        return None
    attrs = tuple(x.strip().split() for x in row[8].rstrip(';').split(';'))  # {gene_id: "AAA", ...}
    attrs = dict((x[0],x[1].strip("\"")) for x in attrs)
    exon_id = attrs.get('exon_id', 'E%d'%Ecounter.next())
    return Exon(id=exon_id,
        gene_id=attrs.get('gene_id',exon_id), gene_name=attrs.get('gene_name',exon_id),
        chrom=row[0], start=int(row[3])-1, end=int(row[4]),
        name=exon_id, score=_score(row[5]), strand=_strand(row[6]),
        transcripts=[attrs.get('transcript_id',exon_id)], exon_number=int(attrs.get('exon_number',1)))

def parse_bed(str line):
    """Parse one BED line. Return False if *line* is empty."""
    cdef list row
    cdef int start,end,strand
    cdef str name
    cdef double score
    if (not line) or (line[0]=='#') or (line[:5]=='track'): return False
    row = line.strip().split('\t')
    chrom = row[0]; start = int(row[1]); end = int(row[2]); name = row[3]
    if len(row) > 4:
        if len(row) > 5: strand = _strand(row[5])
        else: strand = 0
        score = _score(row[4])
    else: score = 0.0
    exon_id = 'E%d'%Ecounter.next()
    return Exon(id=exon_id, gene_id=name, gene_name=name, chrom=chrom, start=start, end=end,
                name=name, score=score, strand=strand, transcripts=[name], exon_number=1)


#########################  Global classes  ##########################


cdef class Counter(object):
    cdef public int n
    def __cinit__(self):
        self.n = 0
    def __call__(self, alignment):
        self.n += 1

cdef class GenomicObject(object):
    cdef public:
        str id,gene_id,gene_name,chrom,name
        int start,end,strand,length,multiplicity
        double score,count,count_rev,rpk,rpk_rev
    def __init__(self,str id='',str gene_id='',str gene_name='',str chrom='',int start=0,int end=0,str name='',
                  double score=0.0,double count=0.0,double count_rev=0.0,double rpk=0.0,double rpk_rev=0.0,
                  int strand=0,int length=0,int multiplicity=1):
        self.id = id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.count = count
        self.count_rev = count_rev
        self.rpk = rpk
        self.rpk_rev = rpk_rev
        self.strand = strand  # 1,-1,0
        self.length = length
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
            name = '|'.join([self.name, other.name]),
            strand = (self.strand + other.strand)/2,
            multiplicity = self.multiplicity + other.multiplicity
        )
    def __repr__(self):
        return "__%s.%s:%d-%d__" % (self.name,self.gene_name,self.start,self.end)

cdef class Exon(GenomicObject):
    cdef public:
        int exon_number
        object transcripts
    cdef:
        double NH,v
        str k
    def __init__(self,int exon_number=0,object transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.exon_number = exon_number
        self.transcripts = transcripts   # list of transcripts it is contained in
        self.length = self.end - self.start
    def __and__(self,other):
        E = GenomicObject.__and__(self,other)
        E.transcripts = set(self.transcripts) | set(other.transcripts)
        return E
    cpdef increment(self,double x,object alignment,bint multiple=False,bint stranded=False):
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
                self.count_rev += x
        else:
            self.count += x

cdef class Transcript(GenomicObject):
    cdef public object exons
    def __init__(self, exons=[], **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons it contains

cdef class Gene(GenomicObject):
    cdef public object exons, transcripts
    def __init__(self, exons=set(),transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons contained
        self.transcripts = transcripts   # list of transcripts contained


#####################  Operations on intervals  #####################


cdef Exon intersect_exons_list(list feats,bint multiple=False):
    """The intersection of a list *feats* of GenomicObjects.
    If *multiple* is True, permits multiplicity: if the same exon E1 is
    given twice, there will be "E1|E1" parts. Otherwise pieces are unique."""
    cdef Exon x,y,f
    if multiple is False:
        feats = list(set(feats))
    if len(feats) == 1:
        #return copy.deepcopy(feats[0])
        # : unavailable in Cython, or must implement the "pickle protocol" with the __reduce__ method
        f = feats[0]
        return Exon(id=f.id,gene_id=f.gene_id,gene_name=f.gene_name,
            chrom=f.chrom,start=f.start,end=f.end,name=f.name,score=f.score,
            strand=f.strand,transcripts=f.transcripts,exon_number=f.exon_number)
    else:
        return reduce(Exon.__and__, feats)

cdef list cobble(list exons,bint multiple=False):
    """Split exons into non-overlapping parts.
    :param multiple: see intersect_exons_list()."""
    cdef list ends, active_exons, cobbled
    cdef Exon e
    cdef tuple a,b
    cdef int i
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


#############################  Counting  ############################


cdef inline double toRPK(double count,double length,double norm_cst):
    return 1000.0 * count / (length * norm_cst)
cdef inline double fromRPK(double rpk,double length,double norm_cst):
    return length * norm_cst * rpk / 1000.0

cdef inline double toSQRPK(double count,double length,double norm_cst):
    return 1000.0 * sqrt(count) / (length * norm_cst)
cdef inline double fromSQRPK(double sqrpk,double length,double norm_cst):
    return pow(length * norm_cst * sqrpk / 1000.0 , 2.0)


cdef int count_reads(list exons,object ckreads,bint multiple,bint stranded) except -1:
    """Adds (#aligned nucleotides/read length) to exon counts.
    Deals with indels, junctions etc.
    :param multiple: divide the count by the NH tag.
    :param standed: for strand-specific protocols, use the strand information."""
    cdef int current_idx, idx2, pos, ali_pos, nexons
    cdef int exon_start, exon_end, shift, op, ali_len, read_len
    cdef object alignment
    cdef Exon E1, E2
    current_idx = 0
    nexons = len(exons)
    for alignment in ckreads:
        if current_idx >= nexons: return 0
        E1 = exons[current_idx]
        exon_end = E1.end
        ali_pos = alignment.pos
        while exon_end <= ali_pos:
            current_idx += 1
            if current_idx >= nexons: return 0
            E1 = exons[current_idx]
            exon_end = E1.end
        idx2 = current_idx
        E2 = exons[idx2]
        exon_start = E2.start
        read_len = alignment.rlen
        ali_len = 0
        for op,shift in alignment.cigar:
            if op in [0,2,3]:  # [BAM_CMATCH,BAM_CDEL,BAM_CREF_SKIP]
                # If read crosses exon left bound
                if ali_pos < exon_start:
                    pos = ali_pos
                    ali_pos = min(exon_start, pos+shift)
                    shift = max(0, pos+shift-exon_start)
                # If read crosses exon right bound, maybe next exon(s)
                while ali_pos+shift >= exon_end:
                    # Score up to exon end, go to exon end, remove from shift and reset ali_len
                    if op == 0:
                        ali_len += exon_end - ali_pos
                        E2.increment(float(ali_len)/float(read_len), alignment, multiple,stranded)
                    shift -= exon_end-ali_pos
                    ali_pos = exon_end
                    ali_len = 0
                    # Next exon
                    idx2 += 1
                    if idx2 >= nexons: return 0
                    E2 = exons[idx2]
                    exon_start = E2.start
                    exon_end = E2.end
                    # If op crosses exon left bound, go to exon start with rest of shift
                    # I op ends before the next exon, go to end of op and reset shift
                    if ali_pos < exon_start:
                        pos = ali_pos
                        ali_pos = min(exon_start, pos+shift)
                        shift = max(0, pos+shift-exon_start)
                # If a bit of op remains overlapping the next exon
                if op == 0:
                    ali_len += shift
                ali_pos += shift   # got to start of next op in prevision for next round
            elif op == 1:  # BAM_CINS
                ali_len += shift;
        E2.increment(float(ali_len)/float(read_len), alignment, multiple,stranded)
    return 0


cdef int get_total_nreads(object sam):
    """Returns the raw total number of alignments (lines) in a *sam* file."""
    cdef str ref
    cdef int length
    cdef Counter Ncounter
    Ncounter = Counter()
    for ref,length in itertools.izip(sam.references,sam.lengths):
        sam.fetch(ref,0,length, callback=Ncounter)
    return Ncounter.n


######################  Expression inference  #######################


cdef inline bint is_in(object feat_class,Exon x,str feat_id):
    """Returns True if Exon *x* is part of the gene/transcript *feat_id*."""
    if feat_class == Transcript:
        return feat_id in x.transcripts
    elif feat_class == Gene:
        return feat_id in x.gene_id.split('|')


cdef list estimate_expression_NNLS(object feat_class,list pieces,list ids,list exons,double norm_cst):
    """Infer gene/transcript expression from exons RPK. Takes Exon instances *pieces*
    and returns for each feature ID in *ids* an instance of *feat_class* with the
    appropriate count and RPK attributes set.

    Builds the exons-transcripts structure matrix:
    Lines are exons, columns are transcripts,
    so that A[i,j]!=0 means 'transcript Tj contains exon Ei'."""
    cdef int n,m,flen,i,j
    cdef double rnorm, fcount, frpk
    cdef str f
    cdef cnp.ndarray[DTYPE_t, ndim=2] A
    cdef cnp.ndarray[DTYPE_t, ndim=1] E, T
    cdef Exon p
    cdef list exs, feats
    n = len(pieces)
    m = len(ids)
    A = zeros((n,m))
    for i,p in enumerate(pieces):
        for j,f in enumerate(ids):
            A[i,j] = 1. if is_in(feat_class,p,f) else 0.
    #--- Build the exons scores vector
    E = asarray([p.rpk for p in pieces])
    #E = asarray([toSQRPK(p.count,p.length,norm_cst) for p in pieces])
    #--- Solve for RPK
    T,rnorm = nnls(A,E)
    #--- Store result in *feat_class* objects
    feats = []
    for i,f in enumerate(ids):
        exs = sorted([e for e in exons if is_in(feat_class,e,f)], key=attrgetter('start','end'))
        flen = sum([p.length for p in pieces if is_in(feat_class,p,f)])
        frpk = T[i]
        fcount = fromRPK(T[i],flen,norm_cst)
        #fSQrpk = T[i]
        #fcount = fromSQRPK(T[i],flen,norm_cst)
        #frpk = toRPK(fcount,flen,norm_cst)
        feats.append(feat_class(name=f, length=flen, rpk=frpk, count=fcount,
                chrom=exs[0].chrom, start=exs[0].start, end=exs[len(exs)-1].end,
                gene_id=exs[0].gene_id, gene_name=exs[0].gene_name))
    return feats


cdef list estimate_expression_raw(object feat_class,list pieces,list ids,list exons,double norm_cst):
    """For each feature ID in *ids*, just sum the score of its components as one
    commonly does for genes from exon counts. Discard ambiguous pieces that are part of
    more than one gene."""
    cdef int flen,i,j
    cdef double fcount, frpk
    cdef str f
    cdef Exon p
    cdef list inner, feats
    feats = []
    for i,f in enumerate(ids):
        exs = sorted([e for e in exons if is_in(feat_class,e,f)], key=attrgetter('start','end'))
        inner = [p for p in pieces if (len(p.gene_id.split('|'))==1 and is_in(feat_class,p,f))]
        flen = sum([p.length for p in inner])
        fcount = sum([p.count for p in inner])
        frpk = toRPK(fcount,flen,norm_cst)
        feats.append(feat_class(name=f, length=flen, rpk=frpk, count=fcount,
                chrom=exs[0].chrom, start=exs[0].start, end=exs[len(exs)-1].end,
                gene_id=exs[0].gene_id, gene_name=exs[0].gene_name))
    return feats


###########################  Main script  ###########################


cdef list partition_chrexons(list chrexons):
    """Partition chrexons in non-overlapping chunks with distinct genes.
    The problem is that exons are sorted wrt start,end, and so the first
    exon of a gene can be separate from the second by exons of other genes
    - from the GTF we don't know how many and how far."""
    cdef int lastend, lastindex, npart, i, lp
    cdef list partition, parts
    cdef dict pinvgenes
    cdef set lastgeneids, toremove
    cdef Exon exon
    cdef str g
    # Cut where disjoint except if the same gene continues
    lastend = chrexons[0].end
    lastgeneids = set(chrexons[0].gene_id)
    lastindex = 0
    partition = []
    pinvgenes = {}  # map {gene_id: partitions it is in}
    npart = 0       # number of partitions
    for i,exon in enumerate(chrexons):
        if (exon.start > lastend) and (exon.gene_id not in lastgeneids):
            lastend = max(exon.end,lastend)
            partition.append((lastindex,i))
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
    npart += 1
    # Merge intervals containing parts of the same gene mixed with others
    toremove = set()
    for g,parts in pinvgenes.iteritems():
        lp = len(parts)
        if lp>1:
            a = parts[0]; b = parts[lp-1]
            partition[b] = (partition[a][0], partition[b][1])
            toremove |= set(xrange(a,b))
    partition = [p for i,p in enumerate(partition) if i not in toremove]
    return partition


def process_chunk(list ckexons,object sam,str chrom,dict options):
    """Distribute counts across transcripts and genes of a chunk *ckexons*
    of non-overlapping exons."""
    cdef int method, lastend
    cdef Exon exon0, g, p
    cdef list pieces, exons, gene_ids, genes, transcripts, el
    cdef tuple es
    cdef dict t2e, e2t, tx_replace
    cdef set filtered
    cdef str t

    #--- Regroup occurrences of the same Exon from a different transcript
    exons = []
    for key,group in itertools.groupby(ckexons, attrgetter('id')):
        # ckexons are sorted by id because chrexons were sorted by chrom,start,end
        exon0 = group.next()
        for g in group:
            exon0.transcripts.append(g.transcripts[0])
        exons.append(exon0)
    gene_ids = list(set(e.gene_id for e in exons))

    #--- Cobble all these intervals
    pieces = cobble(exons)

    #--- Filter out too similar transcripts, e.g. made of the same exons up to 100bp.
    if 1 in options['type']:  # transcripts
        transcript_ids = set()  # full list of remaining transcripts
        t2e = {}                               # map {transcript: [pieces IDs]}
        for p in pieces:
            if p.length < 100: continue        # filter out cobbled pieces of less that read length
            for t in p.transcripts:
                t2e.setdefault(t,[]).append(p.id)
        e2t = {}
        for t,el in t2e.iteritems():
            es = tuple(sorted(el))             # combination of pieces indices
            e2t.setdefault(es,[]).append(t)    # {(pieces IDs combination): [transcripts with same struct]}
        # Replace too similar transcripts by the first of the list, arbitrarily
        tx_replace = dict((badt,tlist[0]) for tlist in e2t.values() for badt in tlist[1:] if len(tlist)>1)
        for p in pieces:
            filtered = set(tx_replace.get(t,t) for t in p.transcripts)
            transcript_ids |= filtered
            p.transcripts = list(filtered)
        transcript_ids = list(transcript_ids)

    #--- Get all reads from this chunk - iterator
    lastend = max(e.end for e in exons)
    ckreads = sam.fetch(chrom, exons[0].start, lastend)

    #--- Count reads in each piece -- from rnacounter.cc
    count_reads(pieces,ckreads,options['multiple'],options['stranded'])

    #--- Calculate RPK
    for p in pieces:
        p.rpk = toRPK(p.count,p.length,options['normalize'])

    #--- Infer gene/transcript counts
    genes = []; transcripts = []
    # Genes - 0
    if 0 in options['type']:
        method = options['method'][0]
        if method == 0:    # raw
            genes = estimate_expression_raw(Gene,pieces,gene_ids,exons,options['normalize'])
        elif method == 1:  # nnls
            genes = estimate_expression_NNLS(Gene,pieces,gene_ids,exons,options['normalize'])
    # Transcripts - 1
    if 1 in options['type']:
        method = options['method'][1]
        if method == 1:    # nnls
            transcripts = estimate_expression_NNLS(Transcript,pieces,transcript_ids,exons,options['normalize'])
        elif method == 0:  # raw
            transcripts = estimate_expression_raw(Transcript,pieces,transcript_ids,exons,options['normalize'])

    #--- Print output
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

    if options['format'] == 'gtf':
        parse = parse_gtf
    elif options['format'] == 'bed':
        parse = parse_bed

    # Initialize
    chrom = ''
    while chrom not in chromosomes:
        exon = None
        while exon is None:
            row = annot.readline().strip()
            exon = parse(row)  # None if not an exon, False if EOF
        if not row: break
        chrom = exon.chrom
    if chrom == '':
        raise ValueError("Reference names in BAM do not correspond to that of the GTF.")
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
                exon = parse(row)  # None if not an exon, False if EOF
            if not row: break
            chrom = exon.chrom
        if len(chrexons) > 0:
            chrexons.sort(key=attrgetter('start','end'))
            partition = partition_chrexons(chrexons)
            # Process chunks
            for (a,b) in partition:
                process_chunk(chrexons[a:b], sam, lastchrom, options)
        lastchrom = chrom

    options['output'].close()
    annot.close()
    sam.close()


########################  Argument parsing  #########################


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
    assert all(x in ["genes","transcripts"] for x in args['--type']), \
        "TYPE must be one of 'genes' or 'transcripts'"
    type_map = {'genes':0, 'transcripts':1, 'exons':2}  # avoid comparing strings later
    args['--type'] = [type_map[x] for x in args['--type']]

    # Same for methods. If given as a list, the length must be that of `--type`,
    # the method at index i will be applied to feature type at index i.
    args['--method'] = [x.lower() for x in args['--method'].split(',')]
    assert len(args['--method']) == len(args['--type'])
    assert all(x in ["raw","nnls","likelihood"] for x in args['--method']), \
        "METHOD must be one of 'raw', 'nnls' or 'likelihood'."
    method_map = {'raw':0, 'nnls':1, 'likelihood':2}  # avoid comparing strings later
    args['--method'] = [method_map[x] for x in args['--method']]
    args['--method'] = dict(zip(args['--type'],args['--method']))

    options = dict((k.lstrip('-').lower(), v) for k,v in args.iteritems())
    return bamname, annotname, options



#----------------------------------------------#
# This code was written by Julien Delafontaine #
# EPFL,BBCF: http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch                       #
#----------------------------------------------#
