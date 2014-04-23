#include "rnacounter.h"

inline int countreads( const samtools::bam1_t *b, void *_d ) {
    size_t i = 0;
    if ((b->core.flag & BAM_FUNMAP) == 0) {
	i = samtools::bam_aux2i(bam_aux_get(b,"NH"));
	if (i==0) i=1;
	(*(std::map< int, size_t >*)_d)[i]++;
    }
    return 0;
}

inline int mapreads( const samtools::bam1_t *b, void *_d ) {
    exon_list *_el = (exon_list*)_d;
    if (_el->current_pos >= _el->size()) return 0;
    int exon_end = (*_el)[_el->current_pos].end,
	ali_pos = b->core.pos;
    while (exon_end <= ali_pos) {
	_el->current_pos++;
	if (_el->current_pos >= _el->size()) return 0;
	exon_end = (*_el)[_el->current_pos].end;
    }
    exon_list::size_type pos2 = _el->current_pos;
    int exon_start = (*_el)[pos2].start,
	read_len = b->core.l_qseq, 
	ali_len = 0;
    for ( uint32_t j = 0; j < b->core.n_cigar; j++ ) {
	int op = bam1_cigar(b)[j]&BAM_CIGAR_MASK,
	    shift = bam1_cigar(b)[j]>>BAM_CIGAR_SHIFT;
	if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP) {
	    if (ali_pos < exon_start) {
		ali_pos = std::min( exon_start, ali_pos+shift );
		shift = std::max( 0, shift-exon_start+ali_pos );
	    }
	    while (ali_pos+shift >= exon_end) {
		if (op == BAM_CMATCH) ali_len += exon_end-ali_pos;
		_el->increment( (double)ali_len/(double)read_len, pos2, bam1_strand(b) );
		ali_pos = exon_end;
		ali_len = 0;
		shift -= exon_end-ali_pos;
		pos2++;
		if (pos2 >= _el->size()) return 0;
		exon_start = (*_el)[pos2].start;
		exon_end = (*_el)[pos2].end;
		if (ali_pos < exon_start) {
		    ali_pos = std::min( exon_start, ali_pos+shift );
		    shift = std::max( 0, shift-exon_start+ali_pos );
		}
	    }
	    ali_pos += shift;
	    if (op == BAM_CMATCH) ali_len += shift;
	}
	if (op == BAM_CINS) ali_len += shift;
    }
    if (ali_len < 1) return 0;
    _el->increment( (double)ali_len/(double)read_len, pos2, bam1_strand(b) );
    return 0;
}



/**************** non-negative least-squares ***************/
#include <libtsnnls/tsnnls.h>
#include <libtsnnls/lsqr.h>

inline void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod ) {
    taucs_ccs_matrix* A = (taucs_ccs_matrix*)prod;
    if ( mode == 0 ) 
	for ( int j = 0; j < A->n; j++) 
	    for ( int i = (A->colptr)[j]; i < (A->colptr[j+1]); i++ )
		y->elements[(A->rowind)[i]] += x->elements[j]*(A->values.d)[i];
    else if ( mode == 1 )
	for ( int j = 0; j < A->n; j++ )
	    for ( int i = 0; i < A->colptr[j+1]-A->colptr[j]; i++ )
		x->elements[j] += y->elements[A->rowind[A->colptr[j]+i]]*A->values.d[A->colptr[j]+i];
    else fprintf(stderr, "Unknown mode: %ld\n", mode );
}

void exon_list::nnls_solve( std::string filename, std::ios_base::openmode mode, int *err ) {
    *err = 0;
/***** rows = #exons, cols = #transcripts *****/
    int rows = size(), cols = txids.size();
    taucs_logfile((char*)"TAUCS.err");
    int nnz = 0;
    std::vector< int > txcounts(cols);
    for ( int n = 0; n < rows; n++ ) {
	nnz += (*this)[n].tmap.size();
	for ( std::vector< int >::iterator I = (*this)[n].tmap.begin();
	      I != (*this)[n].tmap.end();
	      I++ )
	    txcounts[*I]++;
    }
    taucs_ccs_matrix *Avals = taucs_dccs_create(rows, cols, nnz);
    int ncol = 0;
    for ( int n = 0; n < cols; n++ ) {
	Avals->colptr[n] = ncol;
	ncol += txcounts[n];
	txcounts[n] = 0;
    }
    Avals->colptr[cols] = nnz;
    Avals->flags = TAUCS_DOUBLE;
    double norm = 1e3/(double)reads_total;  // normalize to RPKM *= 10^6/rtotal/10^3
    double *bvals_sense, *bvals_anti, *bvals_both; 
    if (stranded) {
	bvals_sense = (double *)(calloc(rows,sizeof(double)));
	bvals_anti = (double *)(calloc(rows,sizeof(double)));
	for ( int n = 0; n < rows; n++ ) {
	    bvals_sense[n] = norm*counts_sense[n];
	    bvals_anti[n]  = norm*counts_antisense[n];
	}
    } else {
	bvals_both = (double *)(calloc(rows,sizeof(double)));
	for ( int n = 0; n < rows; n++ ) bvals_both[n] = norm*counts_both[n];
    }
    for ( int n = 0; n < rows; n++ ) {
	double exlen = (double)((*this)[n].end-(*this)[n].start);
	for ( std::vector< int >::iterator I = (*this)[n].tmap.begin();
	      I != (*this)[n].tmap.end();
	      I++ ) {
	    int pos = txcounts[*I]+Avals->colptr[*I];
	    Avals->rowind[pos] = n;
// should remove 1/2 fragment size (opts) from first and last exon of each transcript 
	    Avals->values.d[pos] = exlen;
	    txcounts[*I]++;
	}
    }
/* DEBUG
   std::cout << "# rows: " << rows << "\n# columns: 1\n";
   for ( int n = 0; n < rows; n++ ) std::cout << bvals[n] << "\n";
   std::cout << "SPARSE\n" << rows << " " << cols << "\n" << nnz << "\n";
   ncol = 1;
   for ( int n = 0; n < nnz; n++ ) {
       if (n == Avals->colptr[ncol]) ncol++;
       std::cout << Avals->rowind[n]+1 << " " << ncol << " " << Avals->values.d[n] << "\n";
   }
    double *xvals, *residual;
    xvals = t_snnls(Avals, bvals, residual, 0.0, 0);
//    xvals = t_lsqr(Avals, bvals);
  if ( lsqr_out->sol_vec->elements == 0 ) {
  std::cerr << "NNLS failed\n";
  *err = 22;
  }
***************************************************/
// ---------- LSQR WRAPPER
    lsqr_input   *lsqr_in;
    lsqr_output  *lsqr_out;
    lsqr_work    *lsqr_work;
    lsqr_func    *lsqr_func;
    alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, rows, cols );
    lsqr_in->num_rows = rows;
    lsqr_in->num_cols = cols;
    lsqr_in->damp_val = 0;
    lsqr_in->rel_mat_err = 0; 
    lsqr_in->rel_rhs_err = 0;
    lsqr_in->cond_lim = 1/sqrt(100*std::numeric_limits<double>::epsilon());
    lsqr_in->max_iter = 4*lsqr_in->num_cols;
    lsqr_in->lsqr_fp_out = NULL;
    if (stranded) {
	for ( int n = 0; n < rows; n++ ) lsqr_in->rhs_vec->elements[n] = bvals_sense[n];
	for ( int n = 0; n < cols; n++ ) lsqr_in->sol_vec->elements[n] = 0;
	lsqr_func->mat_vec_prod = sparse_lsqr_mult;
	lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, Avals );
	print( filename.size() ? filename+std::string(".sense") : std::string(), 
	       lsqr_out->sol_vec->elements, mode );
	free( bvals_sense );
	for ( int n = 0; n < rows; n++ ) lsqr_in->rhs_vec->elements[n] = bvals_anti[n];
	for ( int n = 0; n < cols; n++ ) lsqr_in->sol_vec->elements[n] = 0;
	lsqr_func->mat_vec_prod = sparse_lsqr_mult;
	lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, Avals );
	print( filename.size() ? filename+std::string(".antisense") : std::string(), 
	       lsqr_out->sol_vec->elements, mode );
	free( bvals_anti );
    } else {
	for ( int n = 0; n < rows; n++ ) lsqr_in->rhs_vec->elements[n] = bvals_both[n];
	for ( int n = 0; n < cols; n++ ) lsqr_in->sol_vec->elements[n] = 0;
	lsqr_func->mat_vec_prod = sparse_lsqr_mult;
	lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, Avals );
	print( filename, lsqr_out->sol_vec->elements, mode );
	free( bvals_both );
    }
// -----------------------
    free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
    taucs_dccs_free(Avals);
}

int main( int argc, char **argv )
{
    int err = 0;
    options opts(argc,argv,&err);
    if (err > 0) return err;
    samtools::samfile_t *_bam = samtools::samopen( opts.bamfile.c_str(), "rb", 0 );
    if ( !_bam ) {
	std::cerr << "Could not open " << opts.bamfile << "\n";
	return 10;
    }
    samtools::bam_index_t *_bai = samtools::bam_index_load( opts.bamfile.c_str() );
    if ( !_bai ) {
	std::cerr << "Building index of " << opts.bamfile << "\n";
	samtools::bam_index_build( opts.bamfile.c_str() );
	_bai = samtools::bam_index_load( opts.bamfile.c_str() );
    }
    exon_list exl(&opts);
    if (opts.chroms.empty())
	for ( int cid = 0; cid < _bam->header->n_targets; cid++ )
	    opts.chroms.insert( std::string(_bam->header->target_name[cid]) );
// ---------- count number of reads
    if (opts.normal < 0) {
	std::map< int, size_t > _ntags;
	samtools::bam1_t *b = ((samtools::bam1_t*)calloc(1, sizeof(samtools::bam1_t)));
	while ( samtools::bam_read1( _bam->x.bam, b ) > 0 ) countreads( b, &_ntags );
	exl.update_total( &_ntags );
	free(b);
    }
// ---------- parse bed by chromosome and get corresponding reads
    std::ifstream bedfile( opts.bedfile.c_str() );
    if ( !bedfile.is_open() ) {
	std::cerr << "Could not open " << opts.bedfile << "\n";
	return 11;
    }
    std::ios_base::openmode mode = std::ios_base::out;
    while ( err == 0 && !bedfile.eof() ) {
	std::string oneline;
	std::getline( bedfile, oneline );
	if (exl.append(oneline, opts.chroms)) {
	    int chromid, start, end;
	    samtools::bam_parse_region( _bam->header, exl.chrom.c_str(),
					&chromid, &start, &end );
	    if (chromid < 0) {
		std::cerr << exl.chrom << " not found in " << opts.bamfile << "\n";
		return 12;
	    }
	    exl.reset();
	    samtools::bam_fetch( _bam->x.bam, _bai, chromid, 
				 0, (int)_bam->header->target_len[chromid],
				 &exl, mapreads );
	    exl.nnls_solve( opts.outfile, mode, &err );
	    mode = std::ios_base::app;
	    exl.next();
	}
    }
    samtools::bam_index_destroy( _bai );
    samtools::samclose( _bam );
    return err;
}
