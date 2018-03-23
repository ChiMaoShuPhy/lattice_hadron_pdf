CHROMA=/Users/xiaonuxiong/USQCD/build/chroma
CONFIG=$(CHROMA)/chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I.
LDFLAGS=$(shell $(CONFIG) --ldflags)
LIBS=$(shell $(CONFIG) --libs)

HDRS=kn_meas.h kn_sys.h meson_meas.h nucleon_meas.h spin_matrix.h permutations.h udd.h mom_project.h readwrite.h nucleon.h meson.h myMultiD.h

OBJS=chroma.o \
     kn_meas.o \
     kn_sys.o  \
     meson_meas.o \
     nucleon_meas.o \
     spin_matrix.o \
     permutations.o \
     mom_project.o \
     readwrite.o  \
     meson.o  \
     nucleon.o 


myChroma: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

gammaDriver: gammaDriver.o
	$(CXX) -o $@ $(CXXFLAGS) gammaDriver.o $(LDFLAGS) $(LIBS)

archive:
	tar -zcvf kn.tar.gz *.cc Makefile kn.ini.xml  *.h cfgs/ props/

%.o: %.cc $(HDRS)
	$(CXX) -w $(CXXFLAGS) -c $< 

clean:
	rm -rf $(OBJS) *~

distclean:
	rm -rf *.exe *.o *~ myChroma
