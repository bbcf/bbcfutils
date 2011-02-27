SAMLIBDIR = /archive/epfl/bbcf/bin/samtools
SAMINCLUDE = /archive/epfl/bbcf/bin/
TCLAPINCLUDE = /archive/epfl/bbcf/bin/Peaks/
SQLITELIBDIR = /mnt/common/DevTools/install/Linux/x86_64/sqlite/sqlite-3.7.3/lib 
SQLITEINCLUDE = /mnt/common/DevTools/install/Linux/x86_64/sqlite/sqlite-3.7.3/include/
CXXFLAGS = -O3 -Wall -I$(SAMINCLUDE) -I$(TCLAPINCLUDE) -I$(SQLITEINCLUDE)
LDFLAGS = -L$(SAMLIBDIR) -L$(SQLITELIBDIR) -lsqlite3 -lbam -lz 

SRC = bam2wig.cc bamstat.cc
BIN = bam2wig bamstat

default: $(BIN)

%: %.cc
	$(CXX) -o $* $*.cc $(CXXFLAGS) $(LDFLAGS)

clean:
	$(RM) $(BIN)

cleaner: clean
	$(RM) core *~

