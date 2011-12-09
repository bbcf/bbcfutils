#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
namespace samtools {
#include "samtools/sam.h"
}

/************ 
index  i: nb reads with i hits
index -1: nb of different read starts
index -2: (nb of fwd ali)-(nb of rev ali)
index -i: nb of reads with (i-3) mismatches
*************/
std::map< int, size_t > stats;
int read_length = 0, last_tid = -1, last_start0 = 0, last_start1 = 0;

static inline int samstats( const samtools::bam1_t *b, void *data )
{
    size_t i = 0;
    if (read_length == 0) read_length = b->core.l_qseq;
    if ((b->core.flag & BAM_FUNMAP) == 0) {
	i = samtools::bam_aux2i(bam_aux_get(b,"NH"));
	if (i==0) i=1;
    }
    (*(std::map< int, size_t >*)data)[i]++;
    if ((b->core.flag & BAM_FUNMAP) != 0) return 0;
    i = samtools::bam_aux2i(bam_aux_get(b,"NM"));
    (*(std::map< int, size_t >*)data)[-i-3]++;
    if (b->core.tid != last_tid) {
	last_tid = b->core.tid;
	last_start0 = 0;
	last_start1 = 0;
    }
    int start = b->core.pos+1;
    if (bam1_strand(b)) {
	if (start > last_start1) (*(std::map< int, size_t >*)data)[-1]++;
	last_start1 = start;
	(*(std::map< int, size_t >*)data)[-2]--;
    } else {
	if (start > last_start0) (*(std::map< int, size_t >*)data)[-1]++;
	last_start0 = start;
	(*(std::map< int, size_t >*)data)[-2]++;
    }
    return 0;
}

int main( int argc, char **argv )
{
    samtools::samfile_t *_fs;
    if (argc>1) _fs = samtools::samopen( argv[1], "rb", 0 );
    else 	return 0;
    if ( !_fs ) {
	std::cerr << "Could not open " << argv[1] << "\n";
	return 1;
    }
    samtools::bam_index_t *_in = samtools::bam_index_load( argv[1] );
    if (!_in) {
	std::cerr << "Building index of " << argv[1] << "\n";
	samtools::bam_index_build( argv[1] );
	_in = samtools::bam_index_load( argv[1] );
    }
    int chid = -1, start = 0, end = 0x7fffffff;
    stats.clear();
    long genome_size = 0;
    if (argc>2) {
	samtools::bam_parse_region( _fs->header, argv[2], &chid, &start, &end);
	if (chid < 0) {
	    std::cerr << "Invalid region " << argv[2] << "\n";
	    return 11;
	}
	end = std::min( (int)_fs->header->target_len[chid], end );
	genome_size += (end-start+1);
	samtools::bam_fetch( _fs->x.bam, _in, chid, start, end, &stats, samstats );
    } else {
	for (chid = 0; chid < _fs->header->n_targets; chid++) 
	    genome_size += _fs->header->target_len[chid];
	samtools::bam1_t *b = ((samtools::bam1_t*)calloc(1, sizeof(samtools::bam1_t)));
	while ( samtools::bam_read1( _fs->x.bam, b ) > 0 ) samstats( b, &stats );
    }
    samtools::bam_index_destroy( _in );
    samtools::samclose( _fs );
//-------------- OUTPUT
    std::cout << "Read length " << read_length << "\n";
    std::cout << "Genome size " << genome_size << "\n";
    std::cout << "Nb positions " << stats[-1] << "\n";    
    int ntag = 0, nali = 0;
    std::cout << "Hits Reads\n";
    for ( std::map< int, size_t >::const_iterator I = stats.upper_bound(0); 
	  I != stats.end(); 
	  I++ ) {
	int r = I->second/I->first;
	std::cout << I->first << " " << r << "\n";
	ntag += r;
	nali += I->second;
    }
    std::cout << "Total " << ntag << "\n";
    int nfwd = (nali+stats[-2])/2,
	nrev = (nali-stats[-2])/2;
    std::cout << "Alignments " << nali
	      << " (fwd: " << nfwd << "/rev: " 
	      << nrev << ")\n";
    size_t unmap = 0;
    if ( stats.count(0) ) unmap = stats[0];
    std::cout << "Unmapped " << unmap << "\n";
    std::cout << "Expected coverage " 
	      << std::setiosflags(std::ios::fixed) 
	      << std::setprecision(6)
	      << (double)ntag/(double)(2.0*genome_size)
	      << "\n";
    double act_cov = 0;
    if ( stats.count(-1) ) act_cov = (double)stats[-1]/(double)(2.0*genome_size);
    std::cout << "Actual coverage " << std::setiosflags(std::ios::fixed) 
	      << std::setprecision(6) << act_cov << "\n";
    std::cout << "Mismatches Reads\n";
    for ( std::map< int, size_t >::const_iterator I = stats.upper_bound(-3); 
	  I != stats.begin(); ) {
	I--;
	int r = I->second;
	std::cout << -I->first-3 << " " << r << "\n";
    }
    return 0;
}
