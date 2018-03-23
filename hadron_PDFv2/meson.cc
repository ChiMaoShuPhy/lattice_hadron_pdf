#include "chroma.h"
#include "permutations.h"
#include <iostream>
#include <string>
#include "mom_project.h"
#include "readwrite.h"

#include "meson.h"

using namespace QDP;

int Meson_system::initialize(const multi1d<int>& src,const multi1d<int>& mom, const multi1d<int>& latsize,
                               std::string directory, std::string base)
{
   srce_pt.resize(4);
   lat_size.resize(4);
   p_mom.resize(4);  
    
  /* Set the meson source position */
  srce_pt[0]=src[0];
  srce_pt[1]=src[1];
  srce_pt[2]=src[2];
  srce_pt[3]=src[3];

  /* Set lattice size constants */
  lat_size[0]=latsize[0];
  lat_size[1]=latsize[1];
  lat_size[2]=latsize[2];
  lat_size[3]=latsize[3];
  
      // momentum at the sink
    p_mom[0]=mom[0];
    p_mom[1]=mom[1];
    p_mom[2]=mom[2];
    p_mom[3]=mom[3];
    
    dirname = directory;
    basename = base;
    
    Pxyz.initialize(p_mom,srce_pt);


  return 0;

}

int Meson_system::free()
{
  srce_pt.resize(0);
  lat_size.resize(0);
  p_mom.resize(0);
    
  /* Resize the udd blocks accordingly */
  block.resize(0,0,0,0,0);
  Mblock.resize(0);

  /* Set size of correlators */
  direct.resize( 0 );

  return 0;
  
}

int Meson_system::calcBlock(const LatticePropagator& propQ)
{
  
  quarkPerm = 11;
  
  setSmatrices(propQ);
  
  /* Set size of correlator */
  direct.resize( lat_size[3] );

  /* Resize the meson block accordingly */
  block.resize(4,3,4,3,lat_size[3]);

  if( doChecks ) {
  };
  
  QDPIO::cout << "Meson_system :: Calculating meson block . . . ";
  calcMeson_blocks();
  QDPIO::cout << "Finished." << endl;

  Su.resize(0,0);

  return 0;

}

int Meson_system::calcDirect(const LatticePropagator& propQ)
{
  LatticeComplex temp; 
  /* Set size of correlator */
  direct.resize( lat_size[3] );

  if( doChecks ) {
  };
  
  temp = -trace(adj(propQ) * propQ);
  direct = Pxyz.project_momentum( temp );


  return 0;

}

int Meson_system::calcDirect(const LatticePropagator& propQ, const LatticePropagator& propQbar)
{
  LatticeComplex temp; 
  /* Set size of correlator */
  direct.resize( lat_size[3] );

  if( doChecks ) {
  };
  
  temp = -trace(adj(propQbar) * propQ);
  direct = Pxyz.project_momentum( temp );


  return 0;

}

int Meson_system::calcExchange(const LatticePropagator& propQ,const LatticePropagator& propQbar)
{
  LatticePropagator temp; 


  /* Resize the meson block accordingly */
  Mblock.resize(lat_size[3]);

  if( doChecks ) {
  };
  
  temp = -adj(propQbar) * propQ;
  Mblock = Pxyz.project_momentum(temp);

  return 0;

}

int Meson_system::writeMass()
{
    rw.writeCorrelator(direct,dirname+"meson"+basename,-1, lat_size[3], srce_pt);
    rw.writeAvgCorrelator(direct,dirname+"meson"+basename,-1, lat_size[3], srce_pt[3]);
    return 0;
}

int Meson_system::calcBlock(const LatticePropagator& propQ,const LatticePropagator& propQbar)
{
 
  quarkPerm = 12;
  setSmatrices(propQ,propQbar);

  /* Set size of correlator */
  direct.resize( lat_size[3] );

  /* Resize the meson block accordingly */
  block.resize(4,3,4,3,lat_size[3]);

  if( doChecks ) {
  };

  QDPIO::cout << "Meson_system :: Calculating meson block . . . ";
  calcMeson_blocks();
  QDPIO::cout << "Finished." << endl;

  Su.resize(0,0);
  Ss.resize(0,0);

  return 0;

}

int Meson_system::setSmatrices(const LatticePropagator& propQ)
{

  Su.resize(4,4);
  
  /*  Define all possible color matrices  */
  Su[0][0] = peekSpin(propQ,0,0);
  Su[0][1] = peekSpin(propQ,0,1);
  Su[0][2] = peekSpin(propQ,0,2);
  Su[0][3] = peekSpin(propQ,0,3);
  Su[1][0] = peekSpin(propQ,1,0);
  Su[1][1] = peekSpin(propQ,1,1);
  Su[1][2] = peekSpin(propQ,1,2);
  Su[1][3] = peekSpin(propQ,1,3);
  Su[2][0] = peekSpin(propQ,2,0);
  Su[2][1] = peekSpin(propQ,2,1);
  Su[2][2] = peekSpin(propQ,2,2);
  Su[2][3] = peekSpin(propQ,2,3);
  Su[3][0] = peekSpin(propQ,3,0);
  Su[3][1] = peekSpin(propQ,3,1);
  Su[3][2] = peekSpin(propQ,3,2);
  Su[3][3] = peekSpin(propQ,3,3);

  return 0;

}

int Meson_system::setSmatrices(const LatticePropagator& propQ,const LatticePropagator& propQbar)
{

  Su.resize(4,4);
  Ss.resize(4,4);
  
  /*  Define all possible color matrices  */
  Su[0][0] = peekSpin(propQ,0,0);
  Su[0][1] = peekSpin(propQ,0,1);
  Su[0][2] = peekSpin(propQ,0,2);
  Su[0][3] = peekSpin(propQ,0,3);
  Su[1][0] = peekSpin(propQ,1,0);
  Su[1][1] = peekSpin(propQ,1,1);
  Su[1][2] = peekSpin(propQ,1,2);
  Su[1][3] = peekSpin(propQ,1,3);
  Su[2][0] = peekSpin(propQ,2,0);
  Su[2][1] = peekSpin(propQ,2,1);
  Su[2][2] = peekSpin(propQ,2,2);
  Su[2][3] = peekSpin(propQ,2,3);
  Su[3][0] = peekSpin(propQ,3,0);
  Su[3][1] = peekSpin(propQ,3,1);
  Su[3][2] = peekSpin(propQ,3,2);
  Su[3][3] = peekSpin(propQ,3,3);

  Ss[0][0] = peekSpin(propQbar,0,0);
  Ss[0][1] = peekSpin(propQbar,0,1);
  Ss[0][2] = peekSpin(propQbar,0,2);
  Ss[0][3] = peekSpin(propQbar,0,3);
  Ss[1][0] = peekSpin(propQbar,1,0);
  Ss[1][1] = peekSpin(propQbar,1,1);
  Ss[1][2] = peekSpin(propQbar,1,2);
  Ss[1][3] = peekSpin(propQbar,1,3);
  Ss[2][0] = peekSpin(propQbar,2,0);
  Ss[2][1] = peekSpin(propQbar,2,1);
  Ss[2][2] = peekSpin(propQbar,2,2);
  Ss[2][3] = peekSpin(propQbar,2,3);
  Ss[3][0] = peekSpin(propQbar,3,0);
  Ss[3][1] = peekSpin(propQbar,3,1);
  Ss[3][2] = peekSpin(propQbar,3,2);
  Ss[3][3] = peekSpin(propQbar,3,3);

  return 0;

}

int Meson_system::calcMassFromBlocks()
{
  /*
    Calculates the single Meson masses.
  */
  int a, mu;


  QDPIO::cout << "Meson_system :: Calculating Meson masses . . . ";

  direct = 0.0;
  for(a=0;a<=2;a++){
    for(mu=0;mu<=3;mu++) direct += block[mu][a][mu][a];
    if( doChecks ) {
    }
  }  
  
  /* meson mass correlator */
  rw.writeCorrelator(direct,dirname+"meson"+basename,-1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(direct,dirname+"meson"+basename,-1, lat_size[3], srce_pt[3]);
  
  if( doChecks )
    {
    };

  QDPIO::cout << "Finished." << endl;

  return 0;
}

int Meson_system::calcMeson_blocks()
{
  /* 
     Calculates the meson block
  */
  int alpha, beta, gamma;
  int a, b, c;
  LatticeComplex us;
  
  for (alpha=0;alpha<4;alpha++)
    for (a=0;a<3;a++)
      for (beta=0;beta<4;beta++)
        for (b=0;b<3;b++) {
	  us = 0.0;
	  for (gamma=0;gamma<4;gamma++)
	    for (c=0;c<3;c++){
              switch(quarkPerm){
                  case 11:   
                      us += conj(peekColor(Su[gamma][alpha],c,a)) * peekColor(Su[gamma][beta],c,b);
                      break;
                  case 12:
                      us += conj(peekColor(Ss[gamma][alpha],c,a)) * peekColor(Su[gamma][beta],c,b);
                      break;
              }
	    }
	  block[alpha][a][beta][b]= Pxyz.project_momentum( us );
	}
  
  return 0;
}

