#BASE = /mnt/common/epfl
BASE= /usr
#BASE = /home/leleu/htsstation
PREFIX = $(BASE)/bin
SAMLIBDIR = $(BASE)/lib
SAMINCLUDE = $(BASE)/include
TCLAPINCLUDE = $(BASE)/include
SQLITELIBDIR = /usr/lib64/  
#SQLITELIBDIR = /software/lib64/  
SQLITEINCLUDE = /usr/include/
#SQLITEINCLUDE = /software/include
CXXFLAGS = -O3 -Wall -I$(SAMINCLUDE) -I$(TCLAPINCLUDE) -I$(SQLITEINCLUDE)
LDLIBS = -L$(SAMLIBDIR) -L$(SQLITELIBDIR) -lsqlite3 -lbam -lz -ltsnnls -lsatlas

SRC = bam2wig.cc bamstat.cc rnacounter.cc
BIN = bam2wig bamstat rnacounter S1K S2M SMULTI

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

