CHROMA=/Users/xiaonuxiong/USQCD/build/chroma
CONFIG=$(CHROMA)/chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I. -I/opt/local/include
LDFLAGS=$(shell $(CONFIG) --ldflags) -L/opt/local/lib
LIBS=$(shell $(CONFIG) --libs) -lfftw3_threads -lfftw3

HDRS=

OBJS1= FT_test.o
OBJS2= FT_test_1.o

FT_test: $(OBJS1)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS1) $(LDFLAGS) $(LIBS)

FT_test_1: $(OBJS2)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS2) $(LDFLAGS) $(LIBS)

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf FT_test $(OBJS1) *~
	rm -rf FT_test_1 $(OBJS2) *~
