S1K    = slide_keep.o
S2M    = slide_keep2M.o
SMULTI = slide_multi_keep.o
S2MM    = slide_keep2MM.o


all: S1K S2M SMULTI S2MM

S1K : $(S1K)
	$(CXX) -o S1K $(S1K) $(INCLUDE) $(LIBS)
S2M : $(S2M)
	$(CXX) -o S2M $(S2M) $(INCLUDE) $(LIBS)
S2MM : $(S2MM)
	$(CXX) -o S2MM $(S2MM) $(INCLUDE) $(LIBS)
SMULTI : $(SMULTI)
	$(CXX) -o SMULTI $(SMULTI) $(INCLUDE) $(LIBS)

clean:
	rm -f $(S1K) $(S2M) $(SMULTI) $(S2MM)

$(S1K):    slidetools.h Makefile
$(S2M):    slidetools.h slide_keep2M.cc FixSizeSortedList.h Makefile
$(S2MM):   slidetools.h slide_keep2MM.cc FixSizeSortedList.h Makefile
$(SMULTI): slidetools.h slide_multi_keep.cc Makefile
