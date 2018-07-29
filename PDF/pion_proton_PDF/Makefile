CHROMA=/Users/xiaonuxiong/USQCD/build/chroma
CONFIG=$(CHROMA)/chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I. -I/opt/local/include
LDFLAGS=$(shell $(CONFIG) --ldflags) -L/opt/local/lib
LIBS=$(shell $(CONFIG) --libs) -lfftw3_threads -lfftw3

HDRS=

OBJS1= pion_proton_PDF.o
OBJS2= pion_proton_2_pt.o
OBJS3= pion_proton_2_pt_no_smr.o

## debug
#pion_proton_PDF: $(OBJS1)
#	$(CXX) -g -o $@ $(CXXFLAGS) $(OBJS1) $(LDFLAGS) $(LIBS)

#%.o: %.cc $(HDRS)
#	$(CXX) -g $(CXXFLAGS) -c $<

## release
pion_proton_PDF: $(OBJS1)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS1) $(LDFLAGS) $(LIBS)

pion_proton_2_pt: $(OBJS2)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS2) $(LDFLAGS) $(LIBS)

pion_proton_2_pt_no_smr: $(OBJS3)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS3) $(LDFLAGS) $(LIBS)


%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf pion_proton_PDF $(OBJS1) *~
	rm -rf pion_proton_2_pt $(OBJS2) *~
	rm -rf pion_proton_2_pt_no_smr $(OBJS3) *~