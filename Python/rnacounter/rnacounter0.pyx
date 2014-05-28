
import cython
import pysam
import os, sys, itertools, copy
from numpy import asarray, zeros
from scipy.optimize import nnls

import numpy as np
cimport numpy as cnp
DTYPE = np.double               # fix a datatype for the arrays
ctypedef cnp.double_t DTYPE_t   # assign a corresponding compile-time C type to DTYPE_t


cpdef double _score(x):
    if x == '.': return 0.0
    else: return <double>x
cpdef int _strand(x):
    smap = {'+':1, 1:1, '-':-1, -1:-1, '.':0, 0:0}
    return smap[x]

#Ecounter = itertools.count(1)  # to give unique ids to undefined exons, see parse_gtf()
#def parse_gtf(str row) except -1:
# ...


######################################################################


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
        double score,count,count_rev,rpk
    def __init__(self,str id='',str gene_id='',str gene_name='',str chrom='',int start=0,int end=0,
                  str name='',double score=0.0,double count=0.0,double count_rev=0.0,double rpk=0.0,
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
        return "__%s.%s:%d-%d__" % (self.name,self.gene_name,self.start,self.end)

cdef class Exon(GenomicObject):
    cdef public:
        int exon_number
        object transcripts
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


# Only "def" and not "cdef" because closures (lambdas here) are not supported yet by Cython
def intersect_exons_list(object feats,bint multiple=False):
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
        return reduce(lambda x,y: x&y, feats)

cpdef cobble(object exons,bint multiple=False):
    """Split exons into non-overlapping parts.
    :param multiple: see intersect_exons_list()."""
    cdef object ends, active_exons, cobbled, e, a, b
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


#######################################################################


cpdef partition_chrexons(object chrexons):
    """Partition chrexons in non-overlapping chunks with distinct genes.
    The problem is that exons are sorted wrt start,end, and so the first
    exon of a gene can be separate from the second by exons of other genes
    - from the GTF we don't know how many and how far."""
    cdef int lastend, lastindex, npart, i
    cdef object lastgeneids, pinvgenes, partition, toremove, parts
    cdef Exon exon
    # Cut where disjoint and if the same gene continues
    lastend = chrexons[0].end
    lastgeneids = set([chrexons[0].gene_id])
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
            lastgeneids = set()
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
        if len(parts)>1:
            a = parts[0]; b = parts[-1]
            partition[b] = (partition[a][0], partition[b][1])
            toremove |= set(range(a,b))
    partition = [p for i,p in enumerate(partition) if i not in toremove]

    return partition


cpdef double toRPK(double count,double length,double norm_cst):
    return 1000.0 * count / (length * norm_cst)
cpdef double fromRPK(double rpk,double length,double norm_cst):
    return length * norm_cst * rpk / 1000.


# Only "def" and not "cdef" because closures (lambdas here) are not supported yet by Cython
def estimate_expression(object feat_class,object pieces,object ids,object exons,double norm_cst):
    """Build the exons-transcripts structure matrix:
    Lines are exons, columns are transcripts,
    so that A[i,j]!=0 means 'transcript Tj contains exon Ei'."""
    cdef int n,m,flen,i,j
    cdef double rnorm
    cdef str f
    cdef cnp.ndarray[DTYPE_t, ndim=2] A
    cdef cnp.ndarray[DTYPE_t, ndim=1] E, T
    cdef Exon p
    if feat_class == Gene:
        is_in = lambda Exon x,str g: g in x.gene_id.split('|')
    elif feat_class == Transcript:
        is_in = lambda Exon x,str t: t in x.transcripts
    n = len(pieces)
    m = len(ids)
    A = zeros((n,m))
    for i,p in enumerate(pieces):
        for j,f in enumerate(ids):
            A[i,j] = 1. if is_in(p,f) else 0.
    #--- Build the exons scores vector
    E = asarray([p.rpk for p in pieces])
    #--- Solve for RPK
    T,rnorm = nnls(A,E)
    #--- Store result in *feat_class* objects
    feats = []
    for i,f in enumerate(ids):
        exs = sorted([e for e in exons if is_in(e,f)], key=lambda x:(x.start,x.end))
        flen = sum([p.length for p in pieces if is_in(p,f)])
        feats.append(feat_class(name=f, start=exs[0].start, end=exs[-1].end,
                length=flen, rpk=T[i], count=fromRPK(T[i],flen,norm_cst),
                chrom=exs[0].chrom, gene_id=exs[0].gene_id, gene_name=exs[0].gene_name))
    return feats


cpdef int count_reads(object exons,object ckreads,bint multiple,bint stranded) except -1:
    """Adds (#aligned nucleotides/read length) to exon counts.
    Deals with indels, junctions etc.
    :param multiple: divide the count by the NH tag.
    :param standed: for strand-specific protocols, use the strand information."""
    cdef int current_pos, pos2, ali_pos, nexons
    cdef int exon_start, exon_end, shift, op, ali_len, read_len
    cdef object alignment
    current_pos = 0
    nexons = len(exons)
    for alignment in ckreads:
        if current_pos >= nexons: return 0
        exon_end = exons[current_pos].end
        ali_pos = alignment.pos
        while exon_end <= ali_pos:
            current_pos += 1
            if current_pos >= nexons: return 0
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
                    exons[pos2].increment(float(ali_len)/float(read_len), alignment, multiple,stranded)
                    shift -= exon_end-ali_pos  # remaining part of the read
                    ali_pos = exon_end
                    ali_len = 0
                    pos2 += 1
                    if pos2 >= nexons: return 0
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
        exons[pos2].increment(float(ali_len)/float(read_len), alignment, multiple,stranded)
    return 0


cpdef int get_total_nreads(object sam):
    cdef str ref
    cdef int length
    Ncounter = Counter()
    for ref,length in itertools.izip(sam.references,sam.lengths):
        sam.fetch(ref,0,length, callback=Ncounter)
    return Ncounter.n


#----------------------------------------------#
# This code was written by Julien Delafontaine #
# EPFL,BBCF: http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch                       #
#----------------------------------------------#
