#include "chroma.h"
#include "permutations.h"
#include <iostream>
#include <string>
#include "mom_project.h"
#include "readwrite.h"

#include "kn_sys.h"

using namespace QDP;

int Kaon_Nucleon_system::initialize(const multi1d<int>& srcN, const multi1d<int>& momN,
                                    const multi1d<int>& srcK, const multi1d<int>& momK,
                                    const multi1d<int>& latsize, std::string dir, std::string bas)

{

    dirname = dir;
    basename = bas;
    pK.resize(4);
    pN.resize(4);
    srce_ptK.resize(4);
    srce_ptN.resize(4);
    lat_size.resize(4);
    
  pK[0]=momK[0];pK[1]=momK[1];pK[2]=momK[2];pK[3]=momK[3];
  pN[0]=momN[0];pN[1]=momN[1];pN[2]=momN[2];pN[3]=momN[3];
  srce_ptK[0]=srcK[0];srce_ptK[1]=srcK[1];srce_ptK[2]=srcK[2];srce_ptK[3]=srcK[3];
  srce_ptN[0]=srcN[0];srce_ptN[1]=srcN[1];srce_ptN[2]=srcN[2];srce_ptN[3]=srcN[3];
  lat_size[0]=latsize[0];lat_size[1]=latsize[1];lat_size[2]=latsize[2];lat_size[3]=latsize[3];

  kaon.initialize(srcK,pK,latsize,dirname, basename);
  nucleon.initialize(srcN,pN,latsize,dirname, basename,"Avg");

  return 0;

}

int Kaon_Nucleon_system::free()
{
      nucleon.free();
      kaon.free();

      direct.resize(0);
      exch1.resize(0);
      exch2.resize(0);
      exch3.resize(0);
      
      pK.resize(0);
      pN.resize(0);
      srce_ptK.resize(0);
      srce_ptN.resize(0);
      lat_size.resize(0);
      
      return 0;
}

int Kaon_Nucleon_system::calcExchange(const LatticePropagator& nucleonQ,const LatticePropagator& kaonQ,
                                    const LatticePropagator& kaonS)
{    
    exch1.resize(lat_size[3]);
    exch2.resize(lat_size[3]);
    exch3.resize(lat_size[3]);
    
    QDPIO::cout << "KN_system :: Calculating exchange correlators . . ." << endl;
    nucleon.calcExchange(nucleonQ,kaonQ);
    kaon.calcExchange(nucleonQ,kaonS);
    
    for (int i=0;i<lat_size[3];i++){
        exch1[i] = trace( kaon.Mblock[i] * nucleon.exch1[i] );
        exch2[i] = trace( kaon.Mblock[i] * nucleon.exch2[i] ); 
        exch3[i] = trace( kaon.Mblock[i] * nucleon.exch3[i] ); 
    }
    
    QDPIO::cout << "KN_system :: Finished. . ." << endl;
    return 0;
}

int Kaon_Nucleon_system::calcDirect(const LatticePropagator& nucleonQ,const LatticePropagator& kaonQ,
                                    const LatticePropagator& kaonS)
{
    // This routine calculates the direct part of the KN correlator
    // In doing so, it calculates the kaon and nucleon mass as well
    
    QDPIO::cout << "KN_system :: Calculating direct correlators . . ." << endl;
    
    nucleon.calcDirect(nucleonQ);
    nucleon.writeMass();
    
    kaon.calcDirect(kaonQ,kaonS);
    kaon.writeMass();
    
    direct.resize( lat_size[3] );
    
    direct = kaon.direct * nucleon.direct;  // direct term
    
    QDPIO::cout << "KN_system :: Finished. . ." << endl;
    
    return 0;
}


int Kaon_Nucleon_system::writeCorrelators()
{

  rw.writeCorrelator(direct+exch1,dirname+"kn_I11"+basename,-1, lat_size[3], srce_ptK, srce_ptN);  
  rw.writeAvgCorrelator(direct+exch1,dirname+"kn_I11"+basename,-1, lat_size[3], srce_ptK[3]);
  
  rw.writeCorrelator(direct+exch2+exch3,dirname+"kn_I00"+basename,-1, lat_size[3], srce_ptK, srce_ptN);  // note!!!  I'm off by a minus sign on exch3 !!!
  rw.writeAvgCorrelator(direct+exch2+exch3,dirname+"kn_I00"+basename,-1, lat_size[3], srce_ptK[3]);      //          I don't know why yet!!!!
  
  rw.writeCorrelator(direct+exch2-exch3,dirname+"kn_I10"+basename,-1, lat_size[3], srce_ptK, srce_ptN);  
  rw.writeAvgCorrelator(direct+exch2-exch3,dirname+"kn_I10"+basename,-1, lat_size[3], srce_ptK[3]);

  /*
  rw.writeCorrelator(direct,dirname+"direct"+basename,-1, lat_size[3], srce_ptK, srce_ptN);  
  rw.writeAvgCorrelator(direct,dirname+"direct"+basename,-1, lat_size[3], srce_ptK[3]);
  
  rw.writeCorrelator(exch1,dirname+"exch1"+basename,-1, lat_size[3], srce_ptK, srce_ptN);  
  rw.writeAvgCorrelator(exch1,dirname+"exch1"+basename,-1, lat_size[3], srce_ptK[3]);
  
  rw.writeCorrelator(exch2,dirname+"exch2"+basename,-1, lat_size[3], srce_ptK, srce_ptN);  
  rw.writeAvgCorrelator(exch2,dirname+"exch2"+basename,-1, lat_size[3], srce_ptK[3]);
  
  rw.writeCorrelator(exch3,dirname+"exch3"+basename,-1, lat_size[3], srce_ptK, srce_ptN);  
  rw.writeAvgCorrelator(exch3,dirname+"exch3"+basename,-1, lat_size[3], srce_ptK[3]);
  */
  
  return 0;
}
