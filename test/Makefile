CHROMA=/Users/xiaonuxiong/USQCD/build/chroma
CONFIG=$(CHROMA)/chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I.
LDFLAGS=$(shell $(CONFIG) --ldflags)
LIBS=$(shell $(CONFIG) --libs)

HDRS=

OBJS1= my_test.o
OBJS2= my_test_rgauged.o

my_test: $(OBJS1)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS1) $(LDFLAGS) $(LIBS)

my_test_rgauged: $(OBJS2)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS2) $(LDFLAGS) $(LIBS)

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf my_test $(OBJS1) *~
	rm -rf my_test_rgauged $(OBJS2) *~