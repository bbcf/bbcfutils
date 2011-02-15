/*******************
 g++ -Wall -I/mnt/common/DevTools/install/Linux/x86_64/sqlite/sqlite-3.7.3/include/ -I/home/epfl/bbcf/bin/Peaks/ -I/home/epfl/bbcf/bin/ -O3 -o bam2wig bam2wig.cc -L/home/epfl/bbcf/bin/samtools -L/mnt/common/DevTools/install/Linux/x86_64/sqlite/sqlite-3.7.3/lib -lbam -lz -lsqlite3
******************/
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
    int ntags;
    posh counts;
    posh *scounts;
} samdata;

static const float pseudo_counts(0.5);
static const std::string bed_both("track type=bedGraph visibility=2 name=merged color=10,200,10 windowingFunction=maximum\n");
static const std::string bed_plus("track type=bedGraph visibility=2 name=strand+ color=0,10,200 windowingFunction=maximum\n");
static const std::string bed_minus("track type=bedGraph visibility=2 name=strand- color=200,10,0 windowingFunction=maximum\n");

static struct global_options {
    bool regress, noratio, six, sql;
    int wtpm, mincover, maxhits, maxcnt, cut, cut_ct, chid, start, end, merge;
    std::string ofile, sfile, cfile, chr, chrn;
} opts;

inline double weight_per_tag( const int ntags, 
			      const posh *_cnts=0, const posh *_ctrl=0 ) {
    double weight = 1.0;
//************** 1/[(ntags/10k)/(2*chr size/1Mb)] __average per strand__
    if (opts.wtpm == 0) weight = 2.0e-2*(opts.end-opts.start)/ntags;
    if (opts.wtpm >  0) weight = 1e7/(double)opts.wtpm;
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
    I0--;
    double lc, lastcnt(-1);
    int start = 0, stop = 0;
    sqlite3 *db_fwd, *db_rev, *db_both, *mydb;
    char *sqlErrMsg = 0;
    std::string sql_exec = std::string("create table if not exists ")
	+opts.chrn
	+std::string(" (start integer, end integer, score real)");
    if (opts.merge < 0) {
	std::string fwd = opts.ofile+"fwd";
	if ( sqlite3_open_v2(fwd.c_str(), &db_fwd, 
			     SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL) ) {
	    std::cerr << "Could not open " << fwd << ": " 
		      << sqlite3_errmsg(db_fwd) << "\n";
	    sqlite3_close(db_fwd);
	    return;
	}
	std::string rev = opts.ofile+"rev";
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
	std::string both = opts.ofile+"merged";
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
    sql_exec = std::string("insert  into ")
	+opts.chrn
	+std::string(" (start, end, score) values (?,?,?)");
    if (opts.merge < 0) {
	mydb = db_rev;
	sqlite3_exec( mydb, "begin transaction", NULL, NULL, NULL );
	if ( sqlite3_prepare_v2( mydb, sql_exec.c_str(), 2048, &stmt, &_dummy ) ) {
	    std::cerr << "Prepare error: "  
		      << sqlite3_errmsg(mydb) << "\n";
	    sqlite3_close(mydb);
	    return;
	}
	for ( poshcit I = I0; I != I1; I-- ) {
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
void printbed( posh &counts, const double cntw=1.0 ) {
    if (counts.empty()) return;
    poshcit 
	I0 = counts.lower_bound(-opts.start),
	I1 = counts.lower_bound(-opts.end);
    I0--;
    int lc, lastcnt = -1;
    int start = 0, stop = 0;
    std::ofstream outfile;
    std::ostream* outstr = &std::cout;
    if (opts.ofile.length()) {
	outfile.open( opts.ofile.c_str() );
	if (outfile.is_open()) outstr = &outfile;
	else std::cerr << "Could not open " << opts.ofile 
		       << ", writing to stdout\n";
    }
    if (opts.merge < 0) {
	if (!opts.six) (*outstr) << bed_minus;
	for ( poshcit I = I0; I != I1; I-- ) {
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
    if (!opts.six) {
	if (opts.merge < 0) (*outstr) << bed_plus;
	else                (*outstr) << bed_both;
    }
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
// ****************** NH field = Number of Hits for this tag
    int cur_count = 1;
    if (bam_aux_get(b,"NH")) cur_count = samtools::bam_aux2i(bam_aux_get(b,"NH"));
    if (cur_count <= opts.maxhits || opts.maxhits < 0) {
	float weight = 1.0/cur_count;
	data->ntags++;
// ****************** if cut > 0, set tag size = cut even if longer than read 
	if (opts.cut < 0) opts.cut = b->core.l_qseq;
// ****************** ref is the target set (therefore current set is the control),
// ****************** only record at positions present in ref
	int start = b->core.pos+1, stop = b->core.pos+opts.cut;
	for ( int i = start; i <= stop; i++ ) {
	    if (!data->scounts || data->scounts->count(i)) {
		if (bam1_strand(b)) { // reverse strand
		    if (opts.merge > -1) data->counts[i-opts.merge]+=weight;
		    else                 data->counts[-i]+=weight;
		} else {              // forward strand
		    if (opts.merge > -1) data->counts[i+opts.merge]+=weight;
		    else                 data->counts[i]+=weight;
		}
	    }
	}
    }
    return 1;
}

int main( int argc, char **argv )
{
    try {
	TCLAP::CmdLine cmd( "Reads bam file and creates bedGraph or sqlite with cumulative counts" );
	TCLAP::ValueArg< std::string > of( "o", "outptut", "Output file name (default stdout)", false, "", "string" );
	cmd.add( of ); 
	TCLAP::ValueArg< std::string > sf( "s", "sample", "Sample bam file", true, "", "string" );
	cmd.add( sf ); 
	TCLAP::ValueArg< std::string > cf( "c", "control", "Control bam file", false, "", "string" );
	cmd.add( cf ); 
	TCLAP::ValueArg< std::string > ch( "a", "chromosome", "Chromosome region", true, "", "string" );
	cmd.add( ch ); 
	TCLAP::ValueArg< std::string > cn( "n", "chrname", "Chromosome name", false, "", "string" );
	cmd.add( cn ); 
	TCLAP::ValueArg< int > mincv( "l", "mincover", "Minimum coverage",  false, 0, "int" );
	cmd.add( mincv );
	TCLAP::ValueArg< int > maxct( "t", "maxcounts", "Maximum coverage",  false, 0, "int" );
	cmd.add( maxct );
	TCLAP::ValueArg< int > maxh( "m", "maxhits", "Maximum hit count",  false, 10, "int" );
	cmd.add( maxh );
	TCLAP::ValueArg< int > cut( "q", "cut", "Tags (pseudo-)size",  false, -1, "int" );
	cmd.add( cut );
	TCLAP::ValueArg< int > cut_c( "k", "kut", "Control tags (pseudo-)size",  false, -1, "int" );
	cmd.add( cut_c );
	TCLAP::ValueArg< int > wtpm( "w", "weight", "Normalise by total tag count (per megabase)",  false, -1, "int" );
	cmd.add( wtpm );
	TCLAP::SwitchArg reg( "r", "regress", "Normalize count by regression on control",  false, false );
	cmd.add( reg );
	TCLAP::SwitchArg nrat( "z", "noratio", "Do not compute ratio by control",  false, false );
	cmd.add( nrat );
	TCLAP::ValueArg< int > merge( "p", "merge", "Shift and merge strands",  false, -1, "int" );
	cmd.add( merge );
	TCLAP::SwitchArg sixc( "6", "sixcolumns", "Six columns bed output (default is 4-col bedgraph)",  false, false );
	cmd.add( sixc );
	TCLAP::SwitchArg sql( "d", "sqldb", "Sqlite output",  false, false );
	cmd.add( sql );

	cmd.parse( argc, argv );
	global_options o = {
	    reg.getValue(),
	    nrat.getValue(),
	    sixc.getValue(),
	    sql.getValue(),
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
	if (opts.chrn.empty()) opts.chrn = opts.chr;
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
    } catch( TCLAP::ArgException &e ) {
	std::cerr << "Error: " << e.error() << " " << e.argId() << "\n";
	return 10;
    }
    samtools::samfile_t *_fs = samtools::samopen( opts.sfile.c_str(), "rb", 0 );
    if ( !_fs ) {
	std::cerr << "Could not open " << opts.sfile << "\n";
	return 1;
    }
    samtools::bam_parse_region( _fs->header, opts.chr.c_str(), 
				&opts.chid, &opts.start, &opts.end);
    if (opts.chid < 0) {
	std::cerr << opts.chr << " not found in " << opts.sfile << "\n";
	return 12;
    }
    opts.end = std::min( (int)_fs->header->target_len[opts.chid], opts.end );
    samtools::bam_index_t *_in = samtools::bam_index_load( opts.sfile.c_str() );
    if (!_in) {
	std::cerr << "Building index of " << opts.sfile << "\n";
	samtools::bam_index_build( opts.sfile.c_str() );
	_in = samtools::bam_index_load( opts.sfile.c_str() );
    }
    samdata s_data;
    s_data.scounts = 0; // just making sure...
    samdata c_data;
    c_data.scounts = &s_data.counts;

    samtools::bam_fetch( _fs->x.bam, _in, opts.chid, opts.start, opts.end, 
			 &s_data, accumulate );
    samtools::bam_index_destroy( _in );
    samtools::samclose( _fs );

    double weight = weight_per_tag( s_data.ntags );
    if ( opts.cfile.size() ) {
	_fs = samtools::samopen( opts.cfile.c_str(), "rb", 0 );
	if ( !_fs ) {
	    std::cerr << "Could not open " << opts.cfile << "\n";
	    return 2;
	}
	_in = samtools::bam_index_load( opts.cfile.c_str() );
	if (!_in) {
	    std::cerr << "Building index of " << opts.cfile << "\n";
	    samtools::bam_index_build( opts.cfile.c_str() );
	    _in = samtools::bam_index_load( opts.cfile.c_str() );
	}
	opts.cut = opts.cut_ct;
	samtools::bam_fetch( _fs->x.bam, _in, opts.chid, opts.start, opts.end,
			     &c_data, accumulate );
	samtools::bam_index_destroy( _in );
	samtools::samclose( _fs );
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
    }
    if (opts.sql) createsql( s_data.counts, weight ); 
    else         printbed( s_data.counts, weight ); 
    return 0;
}


