PREFIX = /mnt/common/epfl/bin 
SAMLIBDIR = /mnt/common/epfl/lib
SAMINCLUDE = /mnt/common/epfl/include
TCLAPINCLUDE = /mnt/common/epfl/include
SQLITELIBDIR = /mnt/common/DevTools/install/Linux/x86_64/sqlite/sqlite-3.7.3/lib 
SQLITEINCLUDE = /mnt/common/DevTools/install/Linux/x86_64/sqlite/sqlite-3.7.3/include/
CXXFLAGS = -O3 -Wall -I$(SAMINCLUDE) -I$(TCLAPINCLUDE) -I$(SQLITEINCLUDE)
LOADLIBES = -L$(SAMLIBDIR) -L$(SQLITELIBDIR)
LDLIBS = -lsqlite3 -lbam -lz

SRC = bam2wig.cc bamstat.cc
BIN = bam2wig bamstat S1K S2M SMULTI

all: $(BIN)

S1K:
	cd scanner; $(MAKE) S1K
	mv scanner/S1K .

S2M:
	cd scanner; $(MAKE) S2M
	mv scanner/S2M .

SMULTI:
	cd scanner; $(MAKE) SMULTI
	mv scanner/SMULTI .

install: all
	install $(BIN) $(PREFIX)

clean:
	$(RM) $(BIN)
	cd scanner; $(MAKE) clean

cleaner: clean
	$(RM) core *~

