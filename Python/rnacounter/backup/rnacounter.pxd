
import cython


# Numpy stuff - http://docs.cython.org/src/tutorial/numpy.html
import numpy as np
cimport numpy as cnp
DTYPE = np.double               # fix a datatype for the arrays
ctypedef cnp.double_t DTYPE_t   # assign a corresponding compile-time C type to DTYPE_t


#cpdef inline object parse_gtf(str row)
# -> Error: closures inside cdef functions not yet supported
cpdef inline double _score(double x):
    if str(x) == '.': return 0
    else: return float(x)
cpdef inline int _strand(x):
    smap = {'+':1, 1:1, '-':-1, -1:-1, '.':0, 0:0}
    return smap[x]

cdef class Counter:
    cpdef int n

#cdef class GenomicObject:
#    cpdef int start,end,strand,length,multiplicity
#    cpdef double count,count_rev,rpk,score
#
#cdef class Exon(GenomicObject):
#    #cpdef increment(self,double x,object alignment,bint multiple=*,bint stranded=*):
#    pass
#
#cdef class Transcript(GenomicObject):
#    pass



#cpdef inline intersect_exons_list(feats,bint multiple=*)
# -> Error: closures inside cdef functions not yet supported
#cpdef inline cobble(exons,bint multiple=*)
# -> Error: closures inside cdef functions not yet supported


# process_chrexons()
# -
# process_chunk()
cpdef inline double toRPK(double count, double length, double norm_cst)
cpdef inline double fromRPK(double rpk, double length, double norm_cst)

cpdef inline count_reads(exons,ckreads,bint multiple,bint stranded,double normalize):
    cdef int current_pos,pos2,ali_pos,exon_end,exon_start,read_len,ali_len,op,shift

cpdef inline estimate_epression(feat_class, pieces, ids):
    cdef int n,m,flen
    cdef double rnorm
    cdef cnp.ndarray[DTYPE_t, ndim=2] A
    cdef cnp.ndarray[DTYPE_t, ndim=1] E, T


#rnacounter_main()

