#include "rnacounter.h"

inline int countreads( const samtools::bam1_t *b, void *_d )
{
    size_t i = 0;
    if ((b->core.flag & BAM_FUNMAP) == 0) {
	i = samtools::bam_aux2i(bam_aux_get(b,"NH"));
	if (i==0) i=1;
	(*(std::map< int, size_t >*)_d)[i]++;
    }
    return 0;
}

inline int mapreads( const samtools::bam1_t *b, void *_d ) 
{
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


int exon_list::nnls_solve( std::string filename, 
			   std::ios_base::openmode mode=std::ios_base::out ) {
/***** dim = #exons, cols = #transcripts *****/
    int dim = size(), cols = txids.size();
    taucs_logfile((char*)"TAUCS.err");
    int nnz = 0;
    std::vector< int > txcounts(cols);
    for ( int n = 0; n < dim; n++ ) {
	nnz += (*this)[n].tmap.size();
	for ( std::vector< int >::iterator I = (*this)[n].tmap.begin();
	      I != (*this)[n].tmap.end();
	      I++ )
	    txcounts[*I]++;
    }
    taucs_ccs_matrix *Avals = taucs_dccs_create(dim, cols, nnz);
    int ncol = 0;
    for ( int n = 0; n < cols; n++ ) {
	Avals->colptr[n] = ncol;
	ncol += txcounts[n];
	txcounts[n] = 0;
    }
    Avals->colptr[cols+1] = nnz;
    Avals->flags = TAUCS_DOUBLE;
    double norm = 1e3/(double)reads_total;  // normalize to RPKM *= 10^6/rtotal/10^3
    double *bvals = (double *)(calloc(dim,sizeof(double)));
    for ( int n = 0; n < dim; n++ ) {
	double exlen = (double)((*this)[n].end-(*this)[n].start);
	bvals[n] = norm*counts_antisense[n];
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
    double *xvals, *residual;
    xvals = t_snnls(Avals, bvals, residual, 0.0, 0);
    taucs_ccs_free(Avals);
    if ( xvals > 0 ) {
	print( filename, xvals, mode );
	free(xvals);
	return 0;
    } else {
	std::cerr << "NNLS failed\n";
	return 22;
    }
}

int main( int argc, char **argv )
{
    int err;
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
    exon_list exl;
    exl.stranded = opts.stranded;
// ---------- count number of reads
    std::map< int, size_t > _ntags;
    samtools::bam1_t *b = ((samtools::bam1_t*)calloc(1, sizeof(samtools::bam1_t)));
    while ( samtools::bam_read1( _bam->x.bam, b ) > 0 ) countreads( b, &_ntags );
    exl.update_total( &_ntags );

// ---------- parse bed by chromosome and get corresponding reads
    std::ifstream bedfile( opts.bedfile.c_str() );
    if ( !bedfile.is_open() ) {
	std::cerr << "Could not open " << opts.bedfile << "\n";
	return 11;
    }
    std::string oneline;
    std::getline( bedfile, oneline );
    int chromid, start, end;
    std::ios_base::openmode mode = std::ios_base::out;
    while ( !bedfile.eof() ) {
	if (exl.append(oneline, opts.chrom)) {
	    std::cout << exl.chrom << "\n";
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
	    int err = exl.nnls_solve( opts.outfile, mode );
	    if (err > 0) return err;
	    mode = std::ios_base::app;
	    exl.next();
	}
	std::getline( bedfile, oneline );
    }
    bedfile.close();
    samtools::bam_parse_region( _bam->header, exl.chrom.c_str(), 
				&chromid, &start, &end);
    if (chromid < 0) {
	std::cerr << exl.chrom << " not found in " << opts.bamfile << "\n";
	return 12;
    }
    exl.reset();
    samtools::bam_fetch( _bam->x.bam, _bai, chromid, 
			 0, (int)_bam->header->target_len[chromid],
			 &exl, mapreads );
    exl.nnls_solve( opts.outfile, mode );
    samtools::bam_index_destroy( _bai );
    samtools::samclose( _bam );
    return 0;
}
