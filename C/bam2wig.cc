/***************** STL *************************/
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
/*************** Cmd Line parser ***************/
#include <tclap/CmdLine.h>
/************** samtools: own namespace *******************/
typedef unsigned int uint32_t;
namespace samtools {
#include "samtools/sam.h"
}
/***********************************************/
#include <sqlite3.h> 
/***********************************************/

typedef std::map< int, float > posh;
typedef posh::const_iterator poshcit;
typedef posh::iterator poshit;
typedef struct {
    int strand;
    size_t ntags;
    posh counts;
    posh *scounts;
} samdata;
typedef struct {
    std::string name;
    int start;
    int end;
} chr_region;
std::map< int, chr_region > chr_list;


static const float pseudo_counts(0.5);
static const std::string bed_both("track type=bedGraph visibility=2 name=merged color=10,200,10 windowingFunction=maximum\n");
static const std::string bed_plus("track type=bedGraph visibility=2 name=strand+ color=0,10,200 windowingFunction=maximum\n");
static const std::string bed_minus("track type=bedGraph visibility=2 name=strand- color=200,10,0 windowingFunction=maximum\n");

static struct global_options {
    bool regress, noratio, six, sql, nohds, fragcen;
    long wtpm;
    int mincover, maxhits, maxcnt, cut, cut_ct, chid, start, end, merge;
    std::string ofile, sfile, cfile, chr, chrn;
} opts;

static inline int counttags( const samtools::bam1_t *b, void *data )
{
    size_t i = 0;
    if ((b->core.flag & BAM_FUNMAP) == 0) {
	i = samtools::bam_aux2i(bam_aux_get(b,"NH"));
	if (i==0) i=1;
	(*(std::map< int, size_t >*)data)[i]++;
    }
    return 0;
}

inline double weight_per_tag( const size_t ntags, 
			      const posh *_cnts=0, const posh *_ctrl=0 ) {
    double weight = 1.0;
    if (opts.wtpm == 0L) weight = 1e7/(double)ntags;
    if (opts.wtpm >  0L) weight = 1e7/(double)opts.wtpm;
    if (opts.regress && _ctrl && _ctrl->size()) {
	double scalprod = 0, norm2 = 0;
	for ( poshcit I = _cnts->begin(); I != _cnts->end(); I++ )
	    norm2 += I->second*I->second;
	for ( poshcit I = _ctrl->begin(); I != _ctrl->end(); I++ )
	    scalprod += _cnts->lower_bound(I->first)->second*I->second;
	weight = scalprod/norm2;
    }
    return weight;
}

// ************* sqlite output *************
void createsql( posh &counts, const double cntw=1.0 ) {
    if (counts.empty()) return;
    if (!opts.ofile.length()) {
	std::cerr << "Need an output file name\n";
	return;
    }
    poshcit 
	I0 = counts.lower_bound(-opts.start),
	I1 = counts.lower_bound(-opts.end);
    double lc, lastcnt(-1);
    int start = 0, stop = 0;
    sqlite3 *db_fwd, *db_rev, *db_both, *mydb;
    char *sqlErrMsg = 0;
    std::string sql_exec = std::string("create table if not exists '")+opts.chrn
	+std::string("' (start integer, end integer, score real)");
    if (opts.merge < 0) {
	std::string fwd = opts.ofile+"fwd.sql";
	if ( sqlite3_open_v2(fwd.c_str(), &db_fwd, 
			     SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL) ) {
	    std::cerr << "Could not open " << fwd << ": " 
		      << sqlite3_errmsg(db_fwd) << "\n";
	    sqlite3_close(db_fwd);
	    return;
	}
	std::string rev = opts.ofile+"rev.sql";
	if ( sqlite3_open_v2(rev.c_str(), &db_rev,  
			     SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL) ) {
	    std::cerr << "Could not open " << rev << ": " 
		      << sqlite3_errmsg(db_rev) << "\n";
	    sqlite3_close(db_rev);
	    return;
	}
	if ( sqlite3_exec( db_fwd, sql_exec.c_str(), 
			   NULL, 0, &sqlErrMsg ) ) {
	    std::cerr << "Create table error: " << sqlErrMsg << "\n";
	    sqlite3_free(sqlErrMsg);
	}
	if ( sqlite3_exec( db_rev, sql_exec.c_str(), 
			   NULL, 0, &sqlErrMsg ) ) {
	    std::cerr << "Create table error: " << sqlErrMsg << "\n";
	    sqlite3_free(sqlErrMsg);
	}
    } else {
	std::string both = opts.ofile+"merged.sql";
	if ( sqlite3_open_v2(both.c_str(), &db_both, 
			     SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL) ) {
	    std::cerr << "Could not open " << both << ": " 
		      << sqlite3_errmsg(db_both) << "\n";
	    sqlite3_close(db_both);
	    return;
	}
	if ( sqlite3_exec( db_both, sql_exec.c_str(),
			   NULL, 0, &sqlErrMsg ) ) {
	    std::cerr << "Create table error: " << sqlErrMsg << "\n";
	    sqlite3_free(sqlErrMsg);
	}
    }
    sqlite3_stmt *stmt;
    const char *_dummy;
    sql_exec = std::string("insert into '")+opts.chrn
	+std::string("' (start, end, score) values (?,?,?)");
    if (opts.merge < 0) {
	mydb = db_rev;
	sqlite3_exec( mydb, "begin transaction", NULL, NULL, NULL );
	if ( sqlite3_prepare_v2( mydb, sql_exec.c_str(), 2048, &stmt, &_dummy ) ) {
	    std::cerr << "Prepare error: "  
		      << sqlite3_errmsg(mydb) << "\n";
	    sqlite3_close(mydb);
	    return;
	}
	for ( poshcit I = I0; I != I1; ) {
            I--;
	    lc = I->second*cntw;
	    if (-I->first > stop+1 || std::abs(lastcnt-lc)>1e-6) {
		if (lastcnt >= opts.mincover && 
		    (opts.maxcnt < 1 || lastcnt <= opts.maxcnt)) {
		    sqlite3_bind_int( stmt, 1, start );
		    sqlite3_bind_int( stmt, 2, stop );
		    sqlite3_bind_double( stmt, 3, lastcnt );
		    sqlite3_step( stmt );
		    sqlite3_reset( stmt );
		}
		start = -I->first-1;
		lastcnt = lc;
	    } 
	    stop = -I->first;
	}
	if (lastcnt >= opts.mincover && 
	    (opts.maxcnt < 1 || lastcnt <= opts.maxcnt)) {
	    sqlite3_bind_int( stmt, 1, start );
	    sqlite3_bind_int( stmt, 2, stop );
	    sqlite3_bind_double( stmt, 3, lastcnt );
	    sqlite3_step( stmt );
	}
	sqlite3_finalize( stmt );
	sqlite3_exec( mydb, "commit transaction", NULL, NULL, NULL );
	sqlite3_close( mydb );
    }
    lastcnt = -1;
    start = 0;
    stop = 0;
    if (opts.merge < 0) mydb = db_fwd;
    else mydb = db_both;
    sqlite3_exec( mydb, "begin transaction", NULL, NULL, NULL );
    if ( sqlite3_prepare_v2( mydb, sql_exec.c_str(), 2048, &stmt, &_dummy ) ) {
	std::cerr << "Prepare error: "  
		  << sqlite3_errmsg(mydb) << "\n";
	sqlite3_close(mydb);
	return;
    } 
    I0 = counts.upper_bound(opts.start);
    I1 = counts.upper_bound(opts.end);
    for ( poshcit I = I0; I != I1; I++ ) {
	lc = I->second*cntw;
	if (I->first > stop+1 || std::abs(lastcnt-lc)>1e-6) {
	    if (lastcnt >= opts.mincover && 
		(opts.maxcnt < 1 || lastcnt <= opts.maxcnt)) {
		sqlite3_bind_int( stmt, 1, start );
		sqlite3_bind_int( stmt, 2, stop );
		sqlite3_bind_double( stmt, 3, lastcnt );
		sqlite3_step( stmt );
		sqlite3_reset( stmt );
	    }
	    start = I->first-1;
	    lastcnt = lc;
	} 
	stop = I->first;
    }
    if (lastcnt >= opts.mincover && 
	(opts.maxcnt < 1 || lastcnt <= opts.maxcnt)) {
	    sqlite3_bind_int( stmt, 1, start );
	    sqlite3_bind_int( stmt, 2, stop );
	    sqlite3_bind_double( stmt, 3, lastcnt );
	    sqlite3_step( stmt );
    }
    sqlite3_finalize( stmt );
    sqlite3_exec( mydb, "commit transaction", NULL, NULL, NULL );
    sqlite3_close(mydb);
}

// ************* wiggle track: 0-based starts, 4 cols *************
// *************   bed output: 0-based starts, 6 cols *************
void printbed( posh &counts, const double cntw=1.0, const std::string* header=0,
	       std::ios_base::openmode mode = std::ios_base::out ) {
    if (counts.empty()) return;
    poshcit 
	I0 = counts.lower_bound(-opts.start),
	I1 = counts.lower_bound(-opts.end);
    int lc, lastcnt = -1;
    int start = 0, stop = 0;
    std::ofstream outfile;
    std::ostream* outstr = &std::cout;
    if (opts.ofile.length()) {
	outfile.open( opts.ofile.c_str(), mode );
	if (outfile.is_open()) outstr = &outfile;
	else std::cerr << "Could not open " << opts.ofile 
		       << ", writing to stdout\n";
    }
    if (header) (*outstr) << *header;
    if (opts.merge < 0) {
	for ( poshcit I = I0; I != I1; ) {
            I--;
	    lc = (int)(.5+1e2*I->second*cntw);
	    if (-I->first > stop+1 || lastcnt != lc) {
		if (lastcnt >= 1e2*opts.mincover && 
		    (opts.maxcnt < 1 || lastcnt <= 1e2*opts.maxcnt)) {
		    (*outstr) << opts.chrn << "\t" << start << "\t" << stop << "\t";
		    if (opts.six) (*outstr) << ".\t";
		    (*outstr) << lastcnt*1e-2;
		    if (opts.six) (*outstr) << "\t-";
		    (*outstr) << "\n";
		}
		start = -I->first-1;
		lastcnt = lc;
	    } 
	    stop = -I->first;
	}
	if (lastcnt >= 1e2*opts.mincover && 
	    (opts.maxcnt < 1 || lastcnt <= 1e2*opts.maxcnt)) {
	    (*outstr) << opts.chrn << "\t" << start << "\t" << stop << "\t";
	    if (opts.six) (*outstr) << ".\t";	    
	    (*outstr) << lastcnt*1e-2;
	    if (opts.six) (*outstr) << "\t-";
	    (*outstr) << "\n";
	}
    }
    lastcnt = -1;
    start = 0;
    stop = 0;
    I0 = counts.upper_bound(opts.start);
    I1 = counts.upper_bound(opts.end);
    for ( poshcit I = I0; I != I1; I++ ) {
	lc = (int)(.5+1e2*I->second*cntw);
	if (I->first > stop+1 || lastcnt != lc) {
	    if (lastcnt >= 1e2*opts.mincover && 
		(opts.maxcnt < 1 || lastcnt <= 1e2*opts.maxcnt)) {
		(*outstr) << opts.chrn << "\t" << start << "\t" << stop << "\t";
		if (opts.six) (*outstr) << ".\t";
		(*outstr) << lastcnt*1e-2;
		if (opts.six) (*outstr) << "\t+";
		(*outstr) << "\n";
	    }
	    start = I->first-1;
	    lastcnt = lc;
	} 
	stop = I->first;
    }
    if (lastcnt >= 1e2*opts.mincover && 
	(opts.maxcnt < 1 || lastcnt <= 1e2*opts.maxcnt)) {
	    (*outstr) << opts.chrn << "\t" << start << "\t" << stop << "\t";
	    if (opts.six) (*outstr) << ".\t";	    
	    (*outstr) << lastcnt*1e-2;
	    if (opts.six) (*outstr) << "\t+";
	    (*outstr) << "\n";
    }
    if (outfile.is_open()) outfile.close();
}

inline static int accumulate( const samtools::bam1_t *b, void *d ) {
    samdata *data = (samdata*)d;
    int nh = 1, read_cut = opts.cut, dstr = data->strand;
// ****************** NH field = Number of Hits for this tag
    if (bam_aux_get(b,"NH")) nh = samtools::bam_aux2i(bam_aux_get(b,"NH"));
    if (nh > opts.maxhits && opts.maxhits >= 0) return 1;
    float weight = 1.0/nh;
// ****************** if cut > 0, set tag size = cut even if longer than read 
    if (opts.cut < 0) read_cut = b->core.l_qseq;
    if ((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) return 1;
    if (opts.merge > -1 && (b->core.flag&BAM_FPROPER_PAIR) && !opts.fragcen) {
	read_cut = b->core.isize-opts.merge;
	opts.cut = 1;
	if (read_cut < 1 || bam1_strand(b)) return 1;
    }
//    data->ntags++;
    if (bam1_strand(b) && dstr <= 0) { // ********* reverse strand *********
// ****************** BAM is 0-based, but our array will be 1-based
	int start = b->core.pos+1, stop = start+b->core.l_qseq;
	int st0 = std::max(1,stop-read_cut);
	for ( uint32_t j = 0; j < b->core.n_cigar; j++ ) {
	    int op = bam1_cigar(b)[j]&BAM_CIGAR_MASK,
		shift = bam1_cigar(b)[j]>>BAM_CIGAR_SHIFT;
	    if (op == BAM_CMATCH) {
		int ri0 = std::max(st0,start), ri1 = start+shift;
		if (opts.merge > -1) {
		    ri0 -= opts.merge;
		    ri1 -= opts.merge;
		} else {
		    ri1 = -ri0+1;
		    ri0 = -start-shift+1;
		}
// ****************** ref is the target set (therefore current set is the control),
// ****************** only record at positions present in ref
		for ( int i = ri0; i < ri1; i++ ) 
		    if (!data->scounts || data->scounts->count(i)) 
                        data->counts[i] += weight;
	    }
	    if (op == BAM_CINS || op == BAM_CSOFT_CLIP) 
                st0 += shift;
	    else if (op != BAM_CHARD_CLIP)
                start += shift;
	}
	if (opts.cut > 0) {
	    int ri0 = st0, ri1 = b->core.pos+1;
	    if (opts.merge > -1) {
		ri0 -= opts.merge;
		ri1 -= opts.merge;
	    } else {
		ri0 = -ri1+1;
		ri1 = -st0+1;
	    }
	    for ( int i = ri0; i < ri1; i++ ) 
		if (!data->scounts || data->scounts->count(i)) data->counts[i]+=weight;
	}
    } else if ((!bam1_strand(b)) && dstr >= 0) { // ********* forward strand *********
	int start = b->core.pos+1, stop = start+read_cut;
	if (opts.fragcen) {
	    start = b->core.pos+b->core.isize/2-read_cut/2+1;
	    stop = start+read_cut;
	}
	for ( uint32_t j = 0; j < b->core.n_cigar; j++ ) {
	    int op = bam1_cigar(b)[j]&BAM_CIGAR_MASK,
		shift = bam1_cigar(b)[j]>>BAM_CIGAR_SHIFT;
	    if (op == BAM_CMATCH) {
		int ri0 = start, ri1 = std::min(start+shift, stop);
		if (opts.merge > -1) {
		    ri0 += opts.merge;
		    ri1 += opts.merge;
		}
		for ( int i = ri0; i < ri1; i++ ) 
		    if (!data->scounts || data->scounts->count(i)) 
                        data->counts[i] += weight;
	    }
	    if (op == BAM_CINS || op == BAM_CSOFT_CLIP)
                stop -= shift;
            else if (op != BAM_CHARD_CLIP)
                start += shift;
	    if (op == BAM_CDEL || op == BAM_CREF_SKIP) stop += shift;
	}
	if (opts.cut > 0) {
	    int ri0 = start, ri1 = stop;
	    if (opts.merge > 0) {
		ri0 += opts.merge;
                if ((b->core.flag&BAM_FPROPER_PAIR) == 0) ri1 += opts.merge;
	    }
	    for ( int i = ri0; i < ri1; i++ ) 
		if (!data->scounts || data->scounts->count(i)) data->counts[i]+=weight;
	}
    }
    return 1;
}

int main( int argc, char **argv )
{
    try {
	TCLAP::CmdLine cmd( "Reads bam file and creates bedGraph or sqlite with cumulative counts" );
	TCLAP::ValueArg< std::string > of( "o", "outptut", "Output file prefix (default stdout)", false, "", "string" );
	cmd.add( of ); 
	TCLAP::ValueArg< std::string > sf( "s", "sample", "Sample bam file", true, "", "string" );
	cmd.add( sf ); 
	TCLAP::ValueArg< std::string > cf( "c", "control", "Control bam file", false, "", "string" );
	cmd.add( cf ); 
	TCLAP::ValueArg< std::string > ch( "a", "chromosome", "Chromosome region", false, "", "string" );
	cmd.add( ch ); 
	TCLAP::ValueArg< std::string > cn( "n", "chrname", "Chromosome name (if different from -a)", false, "", "string" );
	cmd.add( cn ); 
	TCLAP::ValueArg< int > mincv( "l", "mincover", "Minimum coverage",  false, 0, "int" );
	cmd.add( mincv );
	TCLAP::ValueArg< int > maxct( "t", "maxcounts", "Maximum coverage",  false, 0, "int" );
	cmd.add( maxct );
	TCLAP::ValueArg< int > maxh( "m", "maxhits", "Maximum number of hits per read",  false, -1, "int" );
	cmd.add( maxh );
	TCLAP::ValueArg< int > cut( "q", "cut", "Tags (pseudo-)size",  false, -1, "int" );
	cmd.add( cut );
	TCLAP::ValueArg< int > cut_c( "k", "kut", "Control tags (pseudo-)size",  false, -1, "int" );
	cmd.add( cut_c );
	TCLAP::ValueArg< long > wtpm( "w", "weight", "If 0: normalise by total tag count*1e-7, if >0: uses 1e-7*w as factor",  false, -1, "int" );
	cmd.add( wtpm );
	TCLAP::SwitchArg reg( "r", "regress", "Normalize count by regression on control",  false );
	cmd.add( reg );
	TCLAP::SwitchArg nrat( "z", "noratio", "Do not compute ratio by control",  false );
	cmd.add( nrat );
	TCLAP::ValueArg< int > merge( "p", "merge", "Shift and merge strands. If >=0 on paired-end, will count whole fragment minus p from each end",  false, -1, "int" );
	cmd.add( merge );
	TCLAP::SwitchArg sixc( "6", "sixcolumns", "Six columns bed output (default is 4-col bedgraph)",  false );
	cmd.add( sixc );
	TCLAP::SwitchArg sql( "d", "sqldb", "Sqlite output",  false );
	cmd.add( sql );
	TCLAP::SwitchArg nohds( "u", "noheaders", "Without bed headers", false );
	cmd.add( nohds );
	TCLAP::SwitchArg fcen( "f", "fragcenter", "For paired-end only: center on fragment midpoint", false );
	cmd.add( fcen );

	cmd.parse( argc, argv );
	global_options o = {
	    reg.getValue(),
	    nrat.getValue(),
	    sixc.getValue(),
	    sql.getValue(),
	    nohds.getValue(),
	    fcen.getValue(),
	    wtpm.getValue(),
	    std::max(0,mincv.getValue()),
	    maxh.getValue(),
	    maxct.getValue(),
	    cut.getValue(),
	    cut_c.getValue(),
	    -1,
	    0,
	    0x7fffffff,
	    merge.getValue(),
	    of.getValue(),
	    sf.getValue(),
	    cf.getValue(),
	    ch.getValue(),
	    cn.getValue()
	};
	opts = o;
	if (opts.cfile.size() 
	    && opts.noratio 
	    && !opts.regress) {
	    std::cerr << "Not using the control, then? (you may want to use -r)\n";
	    opts.cfile.clear();
	}
	if (opts.regress && !opts.cfile.size()) {
	    std::cerr << "Need a control file for regression\n";
	    return 11;
	}
	if (opts.fragcen) opts.merge = 0;
    } catch( TCLAP::ArgException &e ) {
	std::cerr << "Error: " << e.error() << " " << e.argId() << "\n";
	return 10;
    }
/*************** TEST BAM FILE ************/
    samtools::samfile_t *_fs = samtools::samopen( opts.sfile.c_str(), "rb", 0 );
    if ( !_fs ) {
	std::cerr << "Could not open " << opts.sfile << "\n";
	return 1;
    }
    if (opts.chr.empty()) {
	for (int chid = 0; chid < _fs->header->n_targets; chid++) {
	    chr_region cr = { std::string(_fs->header->target_name[chid]),
			      0, _fs->header->target_len[chid] };
	    chr_list[chid] = cr;
	}
    } else {
	samtools::bam_parse_region( _fs->header, opts.chr.c_str(), 
				    &opts.chid, &opts.start, &opts.end);
	if (opts.chid < 0) {
	    std::cerr << opts.chr << " not found in " << opts.sfile << "\n";
	    return 12;
	}
	opts.end = std::min( (int)_fs->header->target_len[opts.chid], opts.end );
	chr_region cr = { std::string(_fs->header->target_name[opts.chid]),
                          opts.start, opts.end };
	chr_list[opts.chid] = cr;
    }
    samtools::bam_index_t *_in = samtools::bam_index_load( opts.sfile.c_str() );
    if (!_in) {
	std::cerr << "Building index of " << opts.sfile << "\n";
	samtools::bam_index_build( opts.sfile.c_str() );
	_in = samtools::bam_index_load( opts.sfile.c_str() );
    }
    samdata s_data;
    s_data.ntags = 0;
    if (opts.wtpm == 0L) {
	std::map< int, size_t > _ntags;
	samtools::bam1_t *b = ((samtools::bam1_t*)calloc(1, sizeof(samtools::bam1_t)));
	while ( samtools::bam_read1( _fs->x.bam, b ) > 0 ) 
	    counttags( b, &_ntags );
	for ( std::map< int, size_t >::const_iterator I = _ntags.begin(); 
	      I != _ntags.end();
	      I++ )
	    s_data.ntags += I->second/I->first;
    }
/*************** CONTROL BAM FILE ************/
    samtools::samfile_t *_cfs = 0;
    samtools::bam_index_t *_cin = 0;
    samdata c_data;
    c_data.ntags = 0;
    if ( opts.cfile.size() ) {
	_cfs = samtools::samopen( opts.cfile.c_str(), "rb", 0 );
	if ( !_cfs ) {
	    std::cerr << "Could not open " << opts.cfile << "\n";
	    return 2;
	}
	_cin = samtools::bam_index_load( opts.cfile.c_str() );
	if (!_cin) {
	    std::cerr << "Building index of " << opts.cfile << "\n";
	    samtools::bam_index_build( opts.cfile.c_str() );
	    _cin = samtools::bam_index_load( opts.cfile.c_str() );
	}
	if (opts.wtpm == 0L) {
	    std::map< int, size_t > _ntags;
	    samtools::bam1_t *b = ((samtools::bam1_t*)calloc(1, sizeof(samtools::bam1_t)));
	    while ( samtools::bam_read1( _cfs->x.bam, b ) > 0 ) 
		counttags( b, &_ntags );
	    for ( std::map< int, size_t >::const_iterator I = _ntags.begin(); 
		  I != _ntags.end();
		  I++ )
		c_data.ntags += I->second/I->first;
	}
    }
    std::ios_base::openmode mode = std::ios_base::out;
    bool chrn_empty = opts.chrn.empty();
    std::vector< int > strset;
    if (opts.merge >= 0 || opts.sql || opts.six) strset.push_back(0);
    else {
        strset.push_back(-1);
        strset.push_back(1);
    }
    for ( std::vector< int >::const_iterator Istr = strset.begin(); Istr != strset.end(); Istr++ ) {
        bool header = !(opts.nohds || opts.sql || opts.six);
        for ( std::map< int, chr_region >::const_iterator I = chr_list.begin();
              I != chr_list.end();
              I++ ) {
            opts.chid = I->first;
            opts.chr = I->second.name;
            if (chrn_empty) opts.chrn = opts.chr;
            opts.start = I->second.start;
            opts.end = I->second.end;
            s_data.counts.clear();
            s_data.scounts = 0; // just making sure...
//            s_data.ntags = 0;
            s_data.strand = *Istr;
            c_data.counts.clear();
            c_data.scounts = &s_data.counts;
//            c_data.ntags = 0;
            c_data.strand = *Istr;

            samtools::bam_fetch( _fs->x.bam, _in, opts.chid, opts.start, opts.end, 
                                 &s_data, accumulate );
	
            double weight = weight_per_tag( s_data.ntags );
            if ( _cfs != 0 ) {
                int old_cut = opts.cut;
                opts.cut = opts.cut_ct;
                samtools::bam_fetch( _cfs->x.bam, _cin, opts.chid, opts.start, opts.end,
                                     &c_data, accumulate );
                weight = weight_per_tag( s_data.ntags, c_data.scounts, &c_data.counts );
                if (!opts.noratio) {
                    double cntw = weight_per_tag(c_data.ntags)/weight;
                    for ( poshit I = s_data.counts.begin(); I != s_data.counts.end(); I++ )
                        I->second += pseudo_counts;
// *********** by construction, keys of control is a subset of keys of counts
                    for ( poshit I = c_data.counts.begin(); I != c_data.counts.end(); I++ )
                        s_data.counts[I->first] /= (pseudo_counts+I->second)*cntw;
                    weight = 1.0;
                }
                opts.cut = old_cut;
            }
            const std::string *_hd = 0;
            if (header) {
                if (*Istr == 0) _hd = &bed_both;
                else if (*Istr > 0) _hd = &bed_plus;
                else if (*Istr < 0) _hd = &bed_minus;
            }
            if (opts.sql) createsql( s_data.counts, weight ); 
            else          printbed( s_data.counts, weight, _hd, mode ); 
            mode = std::ios_base::app;
            header = false;
        }
    }
    samtools::bam_index_destroy( _in );
    samtools::samclose( _fs );
    if ( opts.cfile.size() ) {
        samtools::bam_index_destroy( _cin );
        samtools::samclose( _cfs );
    }
    return 0;
}


