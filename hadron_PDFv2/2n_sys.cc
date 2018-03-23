#include "chroma.h"
#include "permutations.h"
#include <iostream>
#include <string>
#include "mom_project.h"
#include "readwrite.h"

#include "2n_sys.h"

using namespace QDP;

int Two_Nucleon_system::initialize(const multi1d<int>& srcN1, const multi1d<int>& momN1,
                                    const multi1d<int>& srcN2, const multi1d<int>& momN2,
                                    const multi1d<int>& latsize, std::string dir, std::string bas)

{

    dirname = dir;
    basename = bas;
    pN1.resize(4);
    pN2.resize(4);
    srce_ptN1.resize(4);
    srce_ptN2.resize(4);
    lat_size.resize(4);
    
    pN1[0]=momN1[0];pN1[1]=momN1[1];pN1[2]=momN1[2];pN1[3]=momN1[3];
    pN2[0]=momN2[0];pN2[1]=momN2[1];pN2[2]=momN2[2];pN2[3]=momN2[3];
    srce_ptN1[0]=srcN1[0];srce_ptN1[1]=srcN1[1];srce_ptN1[2]=srcN1[2];srce_ptN1[3]=srcN1[3];
    srce_ptN2[0]=srcN2[0];srce_ptN2[1]=srcN2[1];srce_ptN2[2]=srcN2[2];srce_ptN2[3]=srcN2[3];
    lat_size[0]=latsize[0];lat_size[1]=latsize[1];lat_size[2]=latsize[2];lat_size[3]=latsize[3];

    nucleon1.initialize(srcN1,pN1,latsize,dirname, basename,"Avg");
    nucleon2.initialize(srcN2,pN2,latsize,dirname, basename,"Avg");

    return 0;

}

int Two_Nucleon_system::free()
{
      nucleon1.free();
      nucleon2.free();

      direct.resize(0);
      exch1.resize(0);
      exch2.resize(0);
      exch3.resize(0);
      
      pN1.resize(0);
      pN2.resize(0);
      srce_ptN1.resize(0);
      srce_ptN2.resize(0);
      lat_size.resize(0);
      
      return 0;
}

int Two_Nucleon_system::calcExchange(const LatticePropagator& nucleonQ,const LatticePropagator& kaonQ,
                                    const LatticePropagator& kaonS)
{    
    exch1.resize(lat_size[3]);
    exch2.resize(lat_size[3]);
    exch3.resize(lat_size[3]);
    
    QDPIO::cout << "2N_system :: Calculating exchange correlators . . ." << endl;
    nucleon1.calcExchange(nucleonQ,kaonQ);
    nucleon2.calcExchange(nucleonQ,kaonS);
    
    for (int i=0;i<lat_size[3];i++){

    }
    
    QDPIO::cout << "2N_system :: Finished. . ." << endl;
    return 0;
}

int Two_Nucleon_system::calcDirect(const LatticePropagator& nucleonQ,const LatticePropagator& kaonQ,
                                    const LatticePropagator& kaonS)
{
    // This routine calculates the direct part of the 2N correlator
    // In doing so, it calculates the kaon and nucleon mass as well
    
    QDPIO::cout << "2N_system :: Calculating direct correlators . . ." << endl;
    
    nucleon1.calcDirect(nucleonQ);
    nucleon1.writeMass();
    
    nucleon2.calcDirect(kaonQ,kaonS);
    nucleon2.writeMass();
    
    direct.resize( lat_size[3] );
    
//   direct = kaon.direct * nucleon.direct;  // direct term
    
    QDPIO::cout << "2N_system :: Finished. . ." << endl;
    
    return 0;
}


int Two_Nucleon_system::writeCorrelators()
{
  
  return 0;
}
