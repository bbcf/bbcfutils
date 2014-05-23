
import cython


# Numpy stuff - http://docs.cython.org/src/tutorial/numpy.html
import numpy as np
cimport numpy as np
DTYPE = np.double              # fix a datatype for the arrays
ctypedef np.double_t DTYPE_t   # "assign a corresponding compile-time type to DTYPE_t"


#cpdef inline object parse_gtf(str row)
cpdef inline double _score(double x):
    if str(x) == '.': return 0
    else: return float(x)
cpdef inline int _strand(x):
    smap = {'+':1, 1:1, '-':-1, -1:-1, '.':0, 0:0}
    return smap[x]


# Classes...


#cpdef inline intersect_exons_list(feats,bint multiple=*)
#cpdef inline cobble(exons,bint multiple=*)


# process_chrexons()
# -
# process_chunk()
cpdef inline double toRPK(double count, double length, double norm_cst):
    return 1000.0 * count / (length * norm_cst)
cpdef inline double fromRPK(double rpk, double length, double norm_cst):
    return length * norm_cst * rpk / 1000.

cpdef inline count_reads(exons,ckreads,bint multiple,bint stranded,double normalize):
    cdef int current_pos,pos2,ali_pos,exon_end,exon_start,read_len,ali_len,op,shift

cpdef inline estimate_epression(feat_class, pieces, ids):
    cdef int n,m,flen
    cdef double rnorm
    cdef np.ndarray[DTYPE_t, ndim=2] A
    cdef np.ndarray[DTYPE_t, ndim=1] E, T


#rnacounter_main()
cdef class Counter:
    cpdef int n


#cdef class GenomicObject(object):
#    cpdef int start,end,strand,length,multiplicity
#    cpdef double count,count_rev,rpk,score

#    def __cinit__(self, id,gene_id,gene_name,chrom,int start,int end,
#                  name,double score,double count,double count_rev,double rpk,
#                  int strand,int length,seq,int multiplicity):
#        self.id = id
#        self.gene_id = gene_id
#        self.gene_name = gene_name
#        self.chrom = chrom
#        self.start = start
#        self.end = end
#        self.name = name
#        #self.score = score
#        self.count = count
#        self.count_rev = count_rev
#        self.rpk = rpk
#        self.strand = strand
#        self.length = length
#        #self.seq = seq  # sequence
#        self.multiplicity = multiplicity
#    cdef inline object __and__(self, object other)
#    cdef inline str __repr__(self)
