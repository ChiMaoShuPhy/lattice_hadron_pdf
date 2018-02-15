CHROMA=/Users/xiaonuxiong/USQCD/build/chroma
CONFIG=$(CHROMA)/chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I.
LDFLAGS=$(shell $(CONFIG) --ldflags)
LIBS=$(shell $(CONFIG) --libs)

HDRS=

OBJS1= pion_PDF.o
OBJS2= gauge_link.o
OBJS3= multi_shift_test.o

pion_PDF: $(OBJS1)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS1) $(LDFLAGS) $(LIBS)

gauge_link: $(OBJS2)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS3) $(LDFLAGS) $(LIBS)

multi_shift_test: $(OBJS3)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS3) $(LDFLAGS) $(LIBS)

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf pion_PDF $(OBJS1) *~
	rm -rf gauge_link $(OBJS2) *~
	rm -rf multi_shift_test $(OBJS3) *~
