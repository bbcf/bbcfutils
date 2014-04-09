/***************** STL *************************/
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
/*************** Cmd Line parser ***************/
#include <tclap/CmdLine.h>
/************** samtools: own namespace *******************/
typedef unsigned int uint32_t;
namespace samtools {
#include <samtools/sam.h>
}
/**************** non-negative least-squares ***************/
#include <libtsnnls/tsnnls.h>
#undef max
#undef min

typedef struct { int start; int end; std::vector< int > tmap; bool revstrand; } exon;
    
class exon_list : public std::vector< exon > {

public:
// ------ constructor: initialize member vars
    void update_total( std::map< int, size_t >* );
    bool append( std::string, std::string );
    void increment( double, size_type, bool );
    void reset() {
	counts_sense.clear(); counts_antisense.clear(); counts_both.clear();
	if (stranded) {
	    counts_sense.resize(size(),0.0); counts_antisense.resize(size(),0.0);
	} else counts_both.resize(size(),0.0);
	current_pos = 0;
    }
    void next();
    int nnls_solve( std::string, std::ios_base::openmode );
    void print( std::string, double*, std::ios_base::openmode );
    size_type current_pos;
    bool stranded;
    std::string chrom;

private:
    long reads_total;
    std::vector< double > counts_sense;
    std::vector< double > counts_antisense;
    std::vector< double > counts_both;
    std::map< std::string, int > txids;
    std::string next_chrom, next_name;
    exon next_exon;

};

inline void exon_list::update_total( std::map< int, size_t > *_n ) {
    for ( std::map< int, size_t >::const_iterator I = _n->begin(); 
	  I != _n->end();
	  I++ ) reads_total += I->second/I->first;
    _n->clear();
}

inline void exon_list::increment( double x, size_type index, bool read_rev ) {
/* DEBUG   x = 1; */
    if (stranded) {
// ----------- read/exon strand mismatch:
	if (read_rev ^ (*this)[index].revstrand)
	    counts_antisense[index] += x;
	else 
	    counts_sense[index] += x;
    } else counts_both[index] += x;
}

inline bool exon_list::append( std::string line, std::string chr_select ) {
    exon parsed;
    std::string dummy, name, _chr;
    std::stringstream ss(line);
    getline(ss, _chr, '\t');
    if (_chr.compare(chr_select) != 0) return false;
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

inline void exon_list::print( std::string filename, double *x, 
			      std::ios_base::openmode mode=std::ios_base::out ) {
    std::ofstream outfile;
    std::ostream* outstr = &std::cout;
    if (filename.size()) {
	outfile.open( filename.c_str(), mode );
	if (outfile.is_open()) outstr = &outfile;
	else std::cerr << "Could not open " << filename
		       << ", writing to stdout\n";
    }
    for ( std::map< std::string, int >::iterator Im = txids.begin();
	  Im !=  txids.end(); Im++ )
	(*outstr) << Im->first << '\t' << x[Im->second] << "\n";
    if (outfile.is_open()) outfile.close();
}

/*********************** Command-line options ******************/

class options {

public:
    options( int, char**, int* );
    bool stranded;
    std::string bamfile, bedfile, outfile, chrom;

};

options::options( int argc, char **argv, int *err ) {
/********
Add multiple chromosomes (or all), 
write exon-level (bed) output
 ********/
    try {
	TCLAP::CmdLine cmd( "Count reads on transcripts from a genome-level bam file and a non-overlapping set of exons." );
	TCLAP::ValueArg< std::string > of( "o", "outptut", "Output file (default stdout)", false, "", "string" );
	cmd.add( of ); 
	TCLAP::ValueArg< std::string > bf( "b", "bamfile", "Bam file", true, "", "string" );
	cmd.add( bf ); 
	TCLAP::ValueArg< std::string > ef( "e", "exonfile", "Exon bed file", false, "", "string" );
	cmd.add( ef ); 
	TCLAP::ValueArg< std::string > ch( "c", "chromosome", "Chromosome name", false, "", "string" );
	cmd.add( ch );
	TCLAP::SwitchArg st( "s", "stranded", "Compute sense and antisense reads separately",  false );
	cmd.add( st );

	cmd.parse( argc, argv );

	stranded = st.getValue();
	bamfile = bf.getValue();
	bedfile = ef.getValue();
	outfile = of.getValue();
	chrom = ch.getValue();
// DEBUG
	stranded = true;
	bamfile = "/scratch/cluster/monthly/jrougemo/test_rnaQ/Cviwt1_filtered.bam";
	bedfile = "/scratch/cluster/monthly/jrougemo/test_rnaQ/counter/TAIR10_unique_exons.txt";
	outfile = "/scratch/cluster/monthly/jrougemo/test_rnaQ/counter/test_out.txt";
	chrom = "1";
	*err = 0;
    } catch( TCLAP::ArgException &e ) {
	std::cerr << "Error: " << e.error() << " " << e.argId() << "\n";
	*err = 10;
    }
}
