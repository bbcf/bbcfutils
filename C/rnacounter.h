/***************** STL *************************/
#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <limits>
/*************** Cmd Line parser ***************/
#include <tclap/CmdLine.h>
/************** samtools: own namespace *******************/
typedef unsigned int uint32_t;
namespace samtools {
#include <samtools/sam.h>
}
#include <libtsnnls/tsnnls.h>
#include <libtsnnls/lsqr.h>
#undef min
#undef max


#define __FLEN__ 300
/*********************** Command-line options ******************/

class options {

public:
    options( int, char**, int* );
    bool stranded, savecounts;
    int fraglen, normal;
    std::string bamfile, bedfile, outfile;
    std::set< std::string > chroms;

};

options::options( int argc, char **argv, int *err ) {
/********
option : write exon-level (bed) output
 ********/
    try {
	TCLAP::CmdLine cmd( "Count reads on transcripts from a genome-level bam file and a non-overlapping set of exons." );
	TCLAP::ValueArg< std::string > of( "o", "outptut", "Output file (default stdout)", false, "", "string" );
	cmd.add( of ); 
	TCLAP::ValueArg< std::string > bf( "b", "bamfile", "Bam file", true, "", "string" );
	cmd.add( bf ); 
	TCLAP::ValueArg< std::string > ef( "e", "exonfile", "Exon bed file", true, "", "string" );
	cmd.add( ef ); 
	TCLAP::MultiArg< std::string > ch( "c", "chromosome", "Chromosome names", false, "string" );
	cmd.add( ch );
	TCLAP::SwitchArg st( "s", "stranded", "Compute sense and antisense reads separately (false)",  false );
	cmd.add( st );
	TCLAP::SwitchArg si( "i", "intermediate", "Save intermediate values (false)",  false );
	cmd.add( si );
	char buffer [35];
	sprintf(buffer,"Average fragment length (%i)",__FLEN__);
	TCLAP::ValueArg< int > fl( "l", "fraglength", buffer, false, __FLEN__, "int" );
	cmd.add( fl ); 
	TCLAP::ValueArg< int > nm( "n", "normalize", "Normalization constant (total number of reads)", 
				   false, -1, "int" );
	cmd.add( nm ); 

	cmd.parse( argc, argv );

	stranded = st.getValue();
	savecounts = si.getValue();
	fraglen = fl.getValue();
	normal = nm.getValue();
	bamfile = bf.getValue();
	bedfile = ef.getValue();
	outfile = of.getValue();
	std::vector< std::string > _chs = ch.getValue();
	chroms.insert(_chs.begin(), _chs.end());
	*err = 0;
    } catch( TCLAP::ArgException &e ) {
	std::cerr << "Error: " << e.error() << " " << e.argId() << "\n";
	*err = 10;
    }
}

/**************************************************************************************
 exon_list: a list of exons = (start, end, transcipt_set, strand)
            built from parsing a bed-like file:
1       13334   13714   AT1G01030.1     -
1       23145   23415   AT1G01040.1     +
1       23415   24451   AT1G01040.1|AT1G01040.2 +
1       24541   24655   AT1G01040.1|AT1G01040.2 +
 
    The "mapreads" function will be passed to samtools::bam_fetch, call the "increment"
    member function on evry fetched alignment and fill the vectors counts_[sense,antisense,both].
    At the end, "nnls_solve" solves the problem 
             exon_counts = exon_to_transcript_map * transcript_expression
    by least-squares. 
    The exon_to_transcript info (derived from the bedfile's name field as above) is stored in 
    the  
 **************************************************************************************/

typedef struct { int start; int end; std::vector< int > tmap; bool revstrand; } exon;
    
class exon_list : public std::vector< exon > {

public:
    exon_list( const options* );
    void update_total( std::map< int, size_t >* );
    bool append( std::string, std::set< std::string > );
    void increment( double, size_type, bool );
    void reset() {
	counts_sense.clear(); counts_antisense.clear();
	if (stranded) counts_antisense.resize(size(),0.0);
	counts_sense.resize(size(),0.0); 
	current_pos = 0;
    }
    void next();
    void nnls_solve( std::string, std::ios_base::openmode, int* );

    bool stranded, savecounts;
    size_type current_pos;
    int fraglen;
    std::string chrom;

private:
    void print( std::string, double*, std::ios_base::openmode );
    void save_vec_mat( std::string, taucs_ccs_matrix*, double*, double*, std::ios_base::openmode );

    long reads_total;
    std::vector< double > counts_sense; // used for both strands if not stranded
    std::vector< double > counts_antisense; // used only if stranded
    std::map< std::string, int > txids;
    std::string next_chrom, next_name;
    exon next_exon;

};

exon_list::exon_list( const options *o ) {
    stranded = o->stranded;
    savecounts = o->savecounts;
    fraglen = o->fraglen > -1 ? o->fraglen : __FLEN__;
    reads_total = o->normal > -1 ? o->normal : 0;
    chrom = std::string("");
    current_pos = 0;
/*        // DEBUG
	stranded = true;
	bamfile = "/scratch/cluster/monthly/jrougemo/test_rnaQ/Cviwt1_filtered.bam";
	bedfile = "/scratch/cluster/monthly/jrougemo/test_rnaQ/counter/TAIR10_unique_exons_clean.txt";
	outfile = "/scratch/cluster/monthly/jrougemo/test_rnaQ/counter/test_out.txt";
	chroms.insert("1");
*/
}

inline void exon_list::update_total( std::map< int, size_t > *_n ) {
    for ( std::map< int, size_t >::const_iterator I = _n->begin(); 
	  I != _n->end();
	  I++ ) reads_total += I->second/I->first;
    _n->clear();
}

inline void exon_list::increment( double x, size_type index, bool read_rev ) {
// ----------- read/exon strand mismatch:
    if (stranded && (read_rev ^ (*this)[index].revstrand)) counts_antisense[index] += x;
    else                                                   counts_sense[index] += x;
}

inline bool exon_list::append( std::string line, std::set< std::string > chr_select ) {
    if (line.empty() && chr_select.count(chrom) > 0) return true;
    exon parsed;
    std::string dummy, name, _chr;
    std::stringstream ss(line);
    getline(ss, _chr, '\t');    
    if (chr_select.count(_chr) == 0) return false;
    getline(ss, dummy, '\t');
    parsed.start = atoi(dummy.c_str());
    getline(ss, dummy, '\t');
    parsed.end = atoi(dummy.c_str());
    getline(ss, name, '\t');
    getline(ss, dummy, '\t');
    parsed.revstrand = (dummy[0] == '-');
    std::map< std::string, int >::iterator Im;
    if (chrom.empty() || chrom.compare(_chr) == 0) {
	if (chrom.empty()) chrom = _chr;
	std::stringstream ss2(name);
	while (!ss2.eof()) {
	    getline(ss2, dummy, '|');
	    Im = txids.find(dummy);
	    if (Im == txids.end()) 
		Im = txids.insert(Im, std::make_pair(dummy, txids.size()));
	    parsed.tmap.push_back(Im->second);
	}
	push_back( parsed );
	return false;
    }
    next_chrom = _chr;
    next_exon = parsed;
    next_name = name;
    return true;
}

inline void exon_list::next() {
    clear();
    txids.clear();
    chrom = next_chrom;
    std::stringstream ss2(next_name);
    std::string dummy;
    int m = 0;
    while (!ss2.eof()) {
	getline(ss2, dummy, '|');
	txids[dummy] = m;
	next_exon.tmap.push_back(m);
	m++;
    }
    push_back( next_exon );
}



