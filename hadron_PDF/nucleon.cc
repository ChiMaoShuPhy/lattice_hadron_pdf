#include "chroma.h"
#include "permutations.h"
#include <iostream>
#include <string>
#include "mom_project.h"
#include "readwrite.h"

#include "nucleon.h"

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl

using namespace std;
using namespace QDP;

int Nucleon_system::initialize(const multi1d<int>& src,const multi1d<int>& mom, const multi1d<int>& latsize,
                               std::string directory, std::string base, std::string spin)
{
        
    srce_pt.resize(4);
    lat_size.resize(4);
    p_mom.resize(4);
    
    /* Set the source position */
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

    PASS(42);
    
    //    LatticeSpinMatrixD ProjS = zero;
    Complex a;
    std::cout << "#####  " << spin << "!!!!!" << std::endl;
    if (spin=="Avg"){
      std::cout << "here I am" << std::endl;
        a = 0.5;
        pokeSpin(ProjS,a,0,0);
        pokeSpin(ProjS,a,1,1);
    } else if (spin=="Suu") {
        a = 1.0;
        pokeSpin(ProjS,a,0,0);
    } else if (spin=="Sdd") {
        a = 1.0;
        pokeSpin(ProjS,a,1,1);
    } else if (spin=="Sud") {
        a = 1.0;
        pokeSpin(ProjS,a,0,1);
    } else if (spin=="Sdu") {
        a = 1.0;
        pokeSpin(ProjS,a,1,0);
    } else {
      std::cout << "Nothing!!!" << std::endl;
    }
    
    ProjG = 0.0;
    a = 0.7071067811865475244008444;
    pokeSpin(ProjG,a,0,1);
    pokeSpin(ProjG,-a,1,0);

    std::cout << p_mom[0] << " " << p_mom[1] << " " << p_mom[2] << " " << p_mom[3] << std::endl;
    std::cout << srce_pt[0] << " " << srce_pt[1] << " " << srce_pt[2] << " " << srce_pt[3] << std::endl;

    std::cout << "Entering Pxyz initialization" << std::endl;
    Pxyz.initialize(p_mom,srce_pt);
    std::cout << "Pxyz initialization complete!" << std::endl;
    
    return 0;
}

int Nucleon_system::calcExchange(const LatticePropagator& prop1, const LatticePropagator& prop2)
{    
    LatticePropagator Ptemp;
    multi1d<Propagator> P;
    
    exch1.resize(lat_size[3]);
    exch2.resize(lat_size[3]);
    exch3.resize(lat_size[3]);
     
    // why can't I multiply a prop by a multi1d ?????
    
    // note:  I absorb the overall minus sign into these expressions for the exchange terms
    Ptemp = -(ProjG * quarkContract13(prop1,ProjG*prop1)*ProjS + traceSpin(quarkContract13(prop1,ProjG*prop1*ProjG))*ProjS
            +ProjG * transposeSpin(quarkContract34(prop1,ProjS*prop1))*ProjG+ProjS*transposeSpin(quarkContract24(prop1,prop1*ProjG))*ProjG)*prop2;
    exch1 = Pxyz.project_momentum( Ptemp );
    
    Ptemp = -(transposeSpin(quarkContract34(ProjG*prop1*ProjG,ProjS*prop1)) + quarkContract14(ProjS*prop1*ProjG,ProjG*prop1))*prop2;
    exch2 = Pxyz.project_momentum( Ptemp );
    
    Ptemp = -(ProjG * quarkContract14(ProjS*prop1,prop1)*ProjG + ProjG*quarkContract13(ProjG*prop1,prop1)*ProjS
            +quarkContract24(ProjS*prop1,prop1*ProjG)*ProjG+traceSpin(quarkContract13(prop1*ProjG,ProjG*prop1))*ProjS)*prop2;
    exch3 = Pxyz.project_momentum( Ptemp );
               
    return 0;
}

int Nucleon_system::calcDirect(const LatticePropagator& prop1)
{
    LatticeComplex temp;
    
    direct.resize(lat_size[3]);
    
    temp = traceColor(traceSpin(quarkContract13(ProjG*prop1*ProjG,prop1))*traceSpin(ProjS*prop1))
            + trace(quarkContract14(ProjS*prop1*ProjG,ProjG*prop1)*prop1);
    
    direct = Pxyz.project_momentum(temp);
           
    return 0;
}

int Nucleon_system::writeMass()
{
    rw.writeCorrelator(direct,dirname+"nucleon"+basename,-1, lat_size[3], srce_pt);
    rw.writeAvgCorrelator(direct,dirname+"nucleon"+basename,-1, lat_size[3], srce_pt[3]);
    return 0;
}

int Nucleon_system::free()
{
    srce_pt.resize(0);
    lat_size.resize(0);
    p_mom.resize(0);
    
    blockG1pU.resize(0,0,0,0,0,0,0);
    blockG1pD.resize(0,0,0,0,0,0,0);

    corrU.resize( 0 );
    corrD.resize( 0 );
    direct.resize( 0 );
  
    exch1.resize(0);
    exch2.resize(0);
    exch3.resize(0);

    return 0;
  
}

int Nucleon_system::calcBlock(const LatticePropagator& prop)
{
  
  quarkPerm = 111;
  
  setSmatrices(prop);

  /* Resize the udd blocks accordingly */
  blockG1pU.resize(4,3,4,3,4,3,lat_size[3]);
  blockG1pD.resize(4,3,4,3,4,3,lat_size[3]);

  /* Set size of correlators */
  corrU.resize( lat_size[3] );
  corrD.resize( lat_size[3] );

  if( doChecks ) {
  };
  
  QDPIO::cout << "Nucleon_system :: Calculating G1+ nucleon blocks . . . ";
  calcNucleon_blocks(blockG1pU,0,1,0);
  calcNucleon_blocks(blockG1pD,0,1,1);
  QDPIO::cout << "Finished." << endl;

  S1.resize(0,0);
  
  return 0;
}

int Nucleon_system::calcBlock(const LatticePropagator& propU, const LatticePropagator& propD,int qperm)
{
  
  quarkPerm = qperm;
  
  setSmatrices(propU,propD);
  
  /* Resize the udd blocks accordingly */
  blockG1pU.resize(4,3,4,3,4,3,lat_size[3]);
  blockG1pD.resize(4,3,4,3,4,3,lat_size[3]);

  /* Set size of correlators */
  corrU.resize( lat_size[3] );
  corrD.resize( lat_size[3] );

  if( doChecks ) {
  };
  
  QDPIO::cout << "Nucleon_system :: Calculating G1+ nucleon blocks . . . ";
  calcNucleon_blocks(blockG1pU,0,1,0);
  calcNucleon_blocks(blockG1pD,0,1,1);
  QDPIO::cout << "Finished." << endl;

  S1.resize(0,0);
  S2.resize(0,0);

  return 0;
}

int Nucleon_system::calcBlock(const LatticePropagator& propU, const LatticePropagator& propD, 
                               const LatticePropagator& propS)
{
  
  quarkPerm = 123;
  
  setSmatrices(propU,propD,propS);
  
  /* Resize the udd blocks accordingly */
  blockG1pU.resize(4,3,4,3,4,3,lat_size[3]);
  blockG1pD.resize(4,3,4,3,4,3,lat_size[3]);

  /* Set size of correlators */
  corrU.resize( lat_size[3] );
  corrD.resize( lat_size[3] );

  if( doChecks ) {
  };
  
  QDPIO::cout << "Nucleon_system :: Calculating G1+ nucleon blocks . . . ";
  calcNucleon_blocks(blockG1pU,0,1,0);
  calcNucleon_blocks(blockG1pD,0,1,1);
  QDPIO::cout << "Finished." << endl;
  
  S1.resize(0,0);
  S2.resize(0,0);
  S3.resize(0,0);

  return 0;
}

int Nucleon_system::setSmatrices(const LatticePropagator& prop)
{

  S1.resize(4,4);

  S1[0][0] = peekSpin(prop,0,0);
  S1[0][1] = peekSpin(prop,0,1);
  S1[0][2] = peekSpin(prop,0,2);
  S1[0][3] = peekSpin(prop,0,3);
  S1[1][0] = peekSpin(prop,1,0);
  S1[1][1] = peekSpin(prop,1,1);
  S1[1][2] = peekSpin(prop,1,2);
  S1[1][3] = peekSpin(prop,1,3);
  S1[2][0] = peekSpin(prop,2,0);
  S1[2][1] = peekSpin(prop,2,1);
  S1[2][2] = peekSpin(prop,2,2);
  S1[2][3] = peekSpin(prop,2,3);
  S1[3][0] = peekSpin(prop,3,0);
  S1[3][1] = peekSpin(prop,3,1);
  S1[3][2] = peekSpin(prop,3,2);
  S1[3][3] = peekSpin(prop,3,3);
  
  return 0;

}

int Nucleon_system::setSmatrices(const LatticePropagator& prop, const LatticePropagator& propD)
{

  S1.resize(4,4);
  S2.resize(4,4);
    
  /*  Define all possible color matrices  */
  S1[0][0] = peekSpin(prop,0,0);
  S1[0][1] = peekSpin(prop,0,1);
  S1[0][2] = peekSpin(prop,0,2);
  S1[0][3] = peekSpin(prop,0,3);
  S1[1][0] = peekSpin(prop,1,0);
  S1[1][1] = peekSpin(prop,1,1);
  S1[1][2] = peekSpin(prop,1,2);
  S1[1][3] = peekSpin(prop,1,3);
  S1[2][0] = peekSpin(prop,2,0);
  S1[2][1] = peekSpin(prop,2,1);
  S1[2][2] = peekSpin(prop,2,2);
  S1[2][3] = peekSpin(prop,2,3);
  S1[3][0] = peekSpin(prop,3,0);
  S1[3][1] = peekSpin(prop,3,1);
  S1[3][2] = peekSpin(prop,3,2);
  S1[3][3] = peekSpin(prop,3,3);
  
  S2[0][0] = peekSpin(propD,0,0);
  S2[0][1] = peekSpin(propD,0,1);
  S2[0][2] = peekSpin(propD,0,2);
  S2[0][3] = peekSpin(propD,0,3);
  S2[1][0] = peekSpin(propD,1,0);
  S2[1][1] = peekSpin(propD,1,1);
  S2[1][2] = peekSpin(propD,1,2);
  S2[1][3] = peekSpin(propD,1,3);
  S2[2][0] = peekSpin(propD,2,0);
  S2[2][1] = peekSpin(propD,2,1);
  S2[2][2] = peekSpin(propD,2,2);
  S2[2][3] = peekSpin(propD,2,3);
  S2[3][0] = peekSpin(propD,3,0);
  S2[3][1] = peekSpin(propD,3,1);
  S2[3][2] = peekSpin(propD,3,2);
  S2[3][3] = peekSpin(propD,3,3);

  return 0;

}

int Nucleon_system::setSmatrices(const LatticePropagator& prop1, const LatticePropagator& prop2,
                                 const LatticePropagator& prop3)
{

  S1.resize(4,4);
  S2.resize(4,4);
  S3.resize(4,4);
    
  /*  Define all possible color matrices  */
  S1[0][0] = peekSpin(prop1,0,0);
  S1[0][1] = peekSpin(prop1,0,1);
  S1[0][2] = peekSpin(prop1,0,2);
  S1[0][3] = peekSpin(prop1,0,3);
  S1[1][0] = peekSpin(prop1,1,0);
  S1[1][1] = peekSpin(prop1,1,1);
  S1[1][2] = peekSpin(prop1,1,2);
  S1[1][3] = peekSpin(prop1,1,3);
  S1[2][0] = peekSpin(prop1,2,0);
  S1[2][1] = peekSpin(prop1,2,1);
  S1[2][2] = peekSpin(prop1,2,2);
  S1[2][3] = peekSpin(prop1,2,3);
  S1[3][0] = peekSpin(prop1,3,0);
  S1[3][1] = peekSpin(prop1,3,1);
  S1[3][2] = peekSpin(prop1,3,2);
  S1[3][3] = peekSpin(prop1,3,3);
  
  S2[0][0] = peekSpin(prop2,0,0);
  S2[0][1] = peekSpin(prop2,0,1);
  S2[0][2] = peekSpin(prop2,0,2);
  S2[0][3] = peekSpin(prop2,0,3);
  S2[1][0] = peekSpin(prop2,1,0);
  S2[1][1] = peekSpin(prop2,1,1);
  S2[1][2] = peekSpin(prop2,1,2);
  S2[1][3] = peekSpin(prop2,1,3);
  S2[2][0] = peekSpin(prop2,2,0);
  S2[2][1] = peekSpin(prop2,2,1);
  S2[2][2] = peekSpin(prop2,2,2);
  S2[2][3] = peekSpin(prop2,2,3);
  S2[3][0] = peekSpin(prop2,3,0);
  S2[3][1] = peekSpin(prop2,3,1);
  S2[3][2] = peekSpin(prop2,3,2);
  S2[3][3] = peekSpin(prop2,3,3);
  
  S3[0][0] = peekSpin(prop3,0,0);
  S3[0][1] = peekSpin(prop3,0,1);
  S3[0][2] = peekSpin(prop3,0,2);
  S3[0][3] = peekSpin(prop3,0,3);
  S3[1][0] = peekSpin(prop3,1,0);
  S3[1][1] = peekSpin(prop3,1,1);
  S3[1][2] = peekSpin(prop3,1,2);
  S3[1][3] = peekSpin(prop3,1,3);
  S3[2][0] = peekSpin(prop3,2,0);
  S3[2][1] = peekSpin(prop3,2,1);
  S3[2][2] = peekSpin(prop3,2,2);
  S3[2][3] = peekSpin(prop3,2,3);
  S3[3][0] = peekSpin(prop3,3,0);
  S3[3][1] = peekSpin(prop3,3,1);
  S3[3][2] = peekSpin(prop3,3,2);
  S3[3][3] = peekSpin(prop3,3,3);

  return 0;

}

int Nucleon_system::calcMassFromBlocks()
{
  /*
    Calculates the single Nucleon masses using the nucleon blocks
  */
  multi1d<DComplex> invsqrt2( lat_size[3] );
  multi1d<DComplex> eabc( lat_size[3] );
  int a, b, c;

  invsqrt2 = 0.7071067811865475244008444;

  QDPIO::cout << "Nucleon_system :: Calculating Nucleon masses . . . ";

  corrU = 0.0;
  corrD = 0.0;
  for(a=0;a<=2;a++){
    for(b=0;b<=2;b++)
      for(c=0;c<=2;c++){
	if(epsilon(a,b,c) == 0) continue;
	eabc = epsilon(a,b,c);
	corrU += eabc * invsqrt2 * (blockG1pU[0][a][1][b][0][c] - blockG1pU[1][a][0][b][0][c]);
	corrD += eabc * invsqrt2 * (blockG1pD[0][a][1][b][1][c] - blockG1pD[1][a][0][b][1][c]);
	if( doChecks ) {
	}
      };
  };  

   /* blockG1pU */
  rw.writeCorrelator(corrU,dirname+"nucleonG1pU"+basename,-1, lat_size[3], srce_pt);
  rw. writeAvgCorrelator(corrU,dirname+"nucleonG1pU"+basename,-1, lat_size[3], srce_pt[3]);

  /* blockG1pD */
  rw.writeCorrelator(corrD,dirname+"nucleonG1pD"+basename,-1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(corrD,dirname+"nucleonG1pD"+basename,-1, lat_size[3], srce_pt[3]);

  if( doChecks )
    {
      /* blockG1pU */
    };

  QDPIO::cout << "Finished." << endl;

  return 0;
}

multi1d<DComplex> Nucleon_system::scalarMult(double coeff, multi1d<DComplex> md1cm)
{
  multi1d<DComplex> temp(lat_size[3]);
  
  for(int i=0;i<lat_size[3];++i)
    temp[i] = coeff * md1cm[i];

  return temp;
};

int Nucleon_system::calcNucleon_blocks(udd<DComplex>& block, int mu1, int mu2, int mu3)
{
  /* 
     Calculates all the positive parity G1 nucleon blocks.
     Note:  these blocks are open, in that only one contraction over epsilon color indices are performed!!!
  */
  int alpha, beta, gamma;
  int a, b, c;
  double invsqrt2 = 0.7071067811865475244008444;
  
  
  for (alpha=0;alpha<4;alpha++)
    for (a=0;a<3;a++)
      for (beta=0;beta<4;beta++)
	for (b=0;b<3;b++)
	  for(gamma=0;gamma<4;gamma++)
	    for(c=0;c<3;c++)
	      {
		block[alpha][a][beta][b][gamma][c]= scalarMult( invsqrt2 , UDD_permuted(alpha,a,beta,b,gamma,c, mu1,mu2,mu3) 
                                                                         - UDD_permuted(alpha,a,beta,b,gamma,c, mu2,mu1,mu3) );
                if( doChecks )
		  {
		    
		  };
	      }
  
  return 0;
}

multi1d<DComplex> Nucleon_system::UDD_permuted(int alpha, int a, int beta, int b, int gamma, int c, int alphap, int betap, int gammap)
{
  /*
    Returns a ColorMatrix UDD block with zero momentum projection and all {alphap,betap,gammap} permutations summed. 
    Will need to generalize this to arbitrary momentum projection.
  */
  LatticeComplex UDD;
  multi1d<DComplex> temp(lat_size[3]),answer(lat_size[3]);
  int t1,tr1;
  int Nt = lat_size[3];

  UDD = UDD_unpermuted(alpha,a,beta,b,gamma,c, alphap,betap,gammap) + UDD_unpermuted(alpha,a,beta,b,gamma,c, alphap,gammap,betap);

  temp = Pxyz.project_momentum( UDD );
  // now shift the time to have the source time start at zero
  if(srce_pt[3] != 0 ){
    for(t1=0;t1<Nt;t1++)
      if(t1+srce_pt[3] >= Nt){
	tr1=(t1+srce_pt[3]-Nt);
	answer[t1] = -1 * temp[tr1];
      }else{
	tr1=t1+srce_pt[3];
	answer[t1] = temp[tr1];
      }
  } else {
    answer = temp;
  }

  return answer;

}

LatticeComplex Nucleon_system::UDD_unpermuted(int alpha, int a, int beta, int b, int gamma, int c, int alphap, int betap, int gammap){
  /*
    Returns a LatticeColorMatrix UDD block with unpermuted spin indices.
  */

  LatticeComplex SS3;
  int ap;
  int bp;
  int cp;
  int eabcp;

  SS3 = 0;
  for(ap=0;ap<=2;ap++){
    for(bp=0;bp<=2;bp++)
      for(cp=0;cp<=2;cp++){
	eabcp = epsilon(ap,bp,cp);
	if(eabcp == 0) continue;
        switch(quarkPerm){
            case 111:
                SS3 += eabcp * peekColor(S1[alpha][alphap],a,ap) * peekColor(S1[beta][betap],b,bp) * peekColor(S1[gamma][gammap],c,cp);
                break;
            case 112:
                SS3 += eabcp * peekColor(S1[alpha][alphap],a,ap) * peekColor(S1[beta][betap],b,bp) * peekColor(S2[gamma][gammap],c,cp);
                break;
            case 121:
                SS3 += eabcp * peekColor(S1[alpha][alphap],a,ap) * peekColor(S2[beta][betap],b,bp) * peekColor(S1[gamma][gammap],c,cp);
                break;
            case 211:
                SS3 += eabcp * peekColor(S2[alpha][alphap],a,ap) * peekColor(S1[beta][betap],b,bp) * peekColor(S1[gamma][gammap],c,cp);
                break;
            case 123:
                SS3 += eabcp * peekColor(S1[alpha][alphap],a,ap) * peekColor(S2[beta][betap],b,bp) * peekColor(S3[gamma][gammap],c,cp);
                break;
        }    
      }
  };                          /* U                                  D                                D  */

  return SS3;

}

int Nucleon_system::calc_G1pG1p_A1p()
{
  /*
    A1_r1 = 1/sqrt(2) * ( G11 G12  - G12 G11 )
          = sqrt(2) * G11 G12
  */
  multi1d<DComplex> corr( lat_size[3] ), coeff( lat_size[3] );

  QDPIO::cout << "Nucleon_system :: Calculating (G1+G1+)A1+ system . . . ";
  coeff = 2.0 * sqrt(2.0);
  corr = UDD_UDD_permuted(0,1,0, 0,1,1, blockG1pU, blockG1pD); /* G11 G12 Gt11 Gt12 */
  corr = coeff * corr;
  rw.writeCorrelator(corr,dirname+"(G1+G1+)A1+_r1"+basename,1, lat_size[3], srce_pt);
  rw.writeAvgCorrelator(corr,dirname+"(G1+G1+)A1+_r1"+basename,1, lat_size[3], srce_pt[3]);
  
  if( doChecks )
    {

    };
  QDPIO::cout << "Finished." << endl;
  
  return 0;
}

multi1d<DComplex> Nucleon_system::UDD_UDD_permuted(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
					   const udd<DComplex>& UDD1, const udd<DComplex>& UDD2)
{
  multi1d<DComplex> corr( lat_size[3] );

  corr = calcDirect(alpha,beta,gamma,delta,epsilon,phi,UDD1,UDD2);
  
  corr -= calcExchangeU(alpha,beta,gamma,delta,epsilon,phi,UDD1,UDD2);
  corr -= calcExchangeDD1(alpha,beta,gamma,delta,epsilon,phi,UDD1,UDD2);
  corr += calcExchangeDD2(alpha,beta,gamma,delta,epsilon,phi,UDD1,UDD2);
  corr += calcExchangeDD3(alpha,beta,gamma,delta,epsilon,phi,UDD1,UDD2);
  corr -= calcExchangeDD4(alpha,beta,gamma,delta,epsilon,phi,UDD1,UDD2);
  
  return corr;
}

multi1d<DComplex> Nucleon_system::calcDirect(int alpha, int beta, int gamma, int delta, int epsil, int phi, 
					   const udd<DComplex>& UDD1, const udd<DComplex>& UDD2)
{
  multi1d<DComplex> invsqrt2( lat_size[3] );
  multi1d<DComplex> corr11( lat_size[3] );
  multi1d<DComplex> corr22( lat_size[3] );
  multi1d<DComplex> corr12( lat_size[3] );
  multi1d<DComplex> corr21( lat_size[3] );
  multi1d<DComplex> eabc( lat_size[3] );
  int a, b, c;

  invsqrt2 = 0.7071067811865475244008444;

  corr11 = 0.0;
  corr22 = 0.0;
  corr12 = 0.0;
  corr21 = 0.0;
  for(a=0;a<=2;a++){
    for(b=0;b<=2;b++)
      for(c=0;c<=2;c++){
	if(epsilon(a,b,c) == 0) continue;
	eabc = epsilon(a,b,c);
	corr11 += eabc * invsqrt2 * (UDD1[alpha][a][beta][b][gamma][c] - UDD1[beta][a][alpha][b][gamma][c]);
	corr22 += eabc * invsqrt2 * (UDD2[delta][a][epsil][b][phi][c] - UDD2[epsil][a][delta][b][phi][c]);
	corr12 += eabc * invsqrt2 * (UDD2[alpha][a][beta][b][gamma][c] - UDD2[beta][a][alpha][b][gamma][c]);
	corr21 += eabc * invsqrt2 * (UDD1[delta][a][epsil][b][phi][c] - UDD1[epsil][a][delta][b][phi][c]);
      };
  };  
  
  return (corr11 * corr22) - (corr12 * corr21);

}

multi1d<DComplex> Nucleon_system::calcExchangeU(int alphax, int beta, int gamma, int deltax, int epsil, int phi, 
					   const udd<DComplex>& UDD1, const udd<DComplex>& UDD2)
{
  /* 
     Calculates the lone quark exchange:  for example, for neutron, it's the up quark.  For the proton, it's the 
     down quark 
  */
  multi1d<DComplex> corr1( lat_size[3] );
  multi1d<DComplex> corr2( lat_size[3] );
  multi1d<DComplex> inv2( lat_size[3] );
  multi1d<DComplex> eabc( lat_size[3] );
  multi1d<DComplex> eabcp( lat_size[3] );
  int a,ap;
  int b,bp;
  int c,cp;

  corr1 = 0.0;
  corr2 = 0.0;
  inv2 = 0.5;
  
  for(a=0;a<=2;a++)
    for(b=0;b<=2;b++)
      for(c=0;c<=2;c++){
	if(epsilon(a,b,c) == 0) continue;
	eabc = epsilon(a,b,c);
	for(ap=0;ap<=2;ap++){
	  for(bp=0;bp<=2;bp++)
	    for(cp=0;cp<=2;cp++){
	      if(epsilon(ap,bp,cp) == 0) continue;
	      eabcp = epsilon(ap,bp,cp);
	      corr1 += eabc * eabcp * inv2 * ((UDD1[deltax][ap][beta][b][gamma][c]-UDD1[beta][ap][deltax][b][gamma][c]) * (UDD2[alphax][a][epsil][bp][phi][cp] - UDD2[epsil][a][alphax][bp][phi][cp]));
	      corr2 += eabc * eabcp * inv2 * ((UDD1[alphax][ap][epsil][b][phi][c]-UDD1[epsil][ap][alphax][b][phi][c]) * (UDD2[deltax][a][beta][bp][gamma][cp]-UDD2[beta][a][deltax][bp][gamma][cp]));

	    };
	};
      };
  
  return corr1 - corr2;

}

multi1d<DComplex> Nucleon_system::calcExchangeDD1(int alpha, int betax, int gamma, int delta, int epsilx, int phi, 
					   const udd<DComplex>& UDD1, const udd<DComplex>& UDD2)
{
  /* 
     Calculates the double quark exchange:  for example, for neutron, it's one of the down quarks.  For the proton, 
     it's one of the up quarks 
  */
  multi1d<DComplex> corr1( lat_size[3] );
  multi1d<DComplex> corr2( lat_size[3] );
  multi1d<DComplex> inv2( lat_size[3] );
  multi1d<DComplex> eabc( lat_size[3] );
  multi1d<DComplex> eabcp( lat_size[3] );
  int a,ap;
  int b,bp;
  int c,cp;

  corr1 = 0.0;
  corr2 = 0.0;
  inv2 = 0.5;
  
  for(a=0;a<=2;a++)
    for(b=0;b<=2;b++)
      for(c=0;c<=2;c++){
	if(epsilon(a,b,c) == 0) continue;
	eabc = epsilon(a,b,c);
	for(ap=0;ap<=2;ap++){
	  for(bp=0;bp<=2;bp++)
	    for(cp=0;cp<=2;cp++){
	      if(epsilon(ap,bp,cp) == 0) continue;
	      eabcp = epsilon(ap,bp,cp);
	      corr1 += eabc * eabcp * inv2 * ((UDD1[alpha][a][epsilx][bp][gamma][c]-UDD1[epsilx][a][alpha][bp][gamma][c]) * (UDD2[delta][ap][betax][b][phi][cp] - UDD2[betax][ap][delta][b][phi][cp]));
	      corr2 += eabc * eabcp * inv2 * ((UDD1[delta][a][betax][bp][phi][c] - UDD1[betax][a][delta][bp][phi][c]) * (UDD2[alpha][ap][epsilx][b][gamma][cp]-UDD2[epsilx][ap][alpha][b][gamma][cp]));
	    };
	};
      };
  
  return corr1 - corr2;

}

multi1d<DComplex> Nucleon_system::calcExchangeDD2(int alpha, int betax, int gamma, int delta, int epsil, int phix, 
					   const udd<DComplex>& UDD1, const udd<DComplex>& UDD2)
{
  /* 
     Calculates the double quark exchange:  for example, for neutron, it's one of the down quarks.  For the proton, 
     it's one of the up quarks 
  */
  multi1d<DComplex> corr1( lat_size[3] );
  multi1d<DComplex> corr2( lat_size[3] );
  multi1d<DComplex> inv2( lat_size[3] );
  multi1d<DComplex> eabc( lat_size[3] );
  multi1d<DComplex> eabcp( lat_size[3] );
  int a,ap;
  int b,bp;
  int c,cp;

  corr1 = 0.0;
  corr2 = 0.0;
  inv2 = 0.5;
  
  for(a=0;a<=2;a++)
    for(b=0;b<=2;b++)
      for(c=0;c<=2;c++){
	if(epsilon(a,b,c) == 0) continue;
	eabc = epsilon(a,b,c);
	for(ap=0;ap<=2;ap++){
	  for(bp=0;bp<=2;bp++)
	    for(cp=0;cp<=2;cp++){
	      if(epsilon(ap,bp,cp) == 0) continue;
	      eabcp = epsilon(ap,bp,cp);
	      corr1 += eabc * eabcp * inv2 * ((UDD1[alpha][a][phix][bp][gamma][c]-UDD1[phix][a][alpha][bp][gamma][c]) * (UDD2[delta][ap][epsil][bp][betax][c] - UDD2[epsil][ap][delta][bp][betax][c]));
	      corr2 += eabc * eabcp * inv2 * ((UDD1[delta][a][epsil][b][betax][cp] - UDD1[epsil][a][delta][b][betax][cp]) * (UDD2[alpha][ap][phix][b][gamma][cp]-UDD2[phix][ap][alpha][b][gamma][cp]));
	    };
	};
      };
  
  return corr1-corr2;

}


multi1d<DComplex> Nucleon_system::calcExchangeDD3(int alpha, int beta, int gammax, int delta, int epsilx, int phi, 
					   const udd<DComplex>& UDD1, const udd<DComplex>& UDD2)
{
  /* 
     Calculates the double quark exchange:  for example, for neutron, it's one of the down quarks.  For the proton, 
     it's one of the up quarks 
  */
  multi1d<DComplex> corr1( lat_size[3] );
  multi1d<DComplex> corr2( lat_size[3] );
  multi1d<DComplex> inv2( lat_size[3] );
  multi1d<DComplex> eabc( lat_size[3] );
  multi1d<DComplex> eabcp( lat_size[3] );
  int a,ap;
  int b,bp;
  int c,cp;

  corr1 = 0.0;
  corr2 = 0.0;
  inv2 = 0.5;
  
  for(a=0;a<=2;a++)
    for(b=0;b<=2;b++)
      for(c=0;c<=2;c++){
	if(epsilon(a,b,c) == 0) continue;
	eabc = epsilon(a,b,c);
	for(ap=0;ap<=2;ap++){
	  for(bp=0;bp<=2;bp++)
	    for(cp=0;cp<=2;cp++){
	      if(epsilon(ap,bp,cp) == 0) continue;
	      eabcp = epsilon(ap,bp,cp);
	      corr1 += eabc * eabcp * inv2 * ((UDD1[alpha][a][beta][b][epsilx][cp]-UDD1[beta][a][alpha][b][epsilx][cp]) * (UDD2[delta][ap][gammax][b][phi][cp] - UDD2[gammax][ap][delta][b][phi][cp]));
	      corr2 += eabc * eabcp * inv2 * ((UDD1[delta][a][gammax][bp][phi][c] - UDD1[gammax][a][delta][bp][phi][c]) * (UDD2[alpha][ap][beta][bp][epsilx][c]-UDD2[beta][ap][alpha][bp][epsilx][c]));
	    };
	};
      };
  
  return corr1 - corr2;

}

multi1d<DComplex> Nucleon_system::calcExchangeDD4(int alpha, int beta, int gammax, int delta, int epsil, int phix, 
					   const udd<DComplex>& UDD1, const udd<DComplex>& UDD2)
{
  /* 
     Calculates the double quark exchange:  for example, for neutron, it's one of the down quarks.  For the proton, 
     it's one of the up quarks 
  */
  multi1d<DComplex> corr1( lat_size[3] );
  multi1d<DComplex> corr2( lat_size[3] );
  multi1d<DComplex> inv2( lat_size[3] );
  multi1d<DComplex> eabc( lat_size[3] );
  multi1d<DComplex> eabcp( lat_size[3] );
  int a,ap;
  int b,bp;
  int c,cp;

  corr1 = 0.0;
  corr2 = 0.0;
  inv2 = 0.5;
  
  for(a=0;a<=2;a++)
    for(b=0;b<=2;b++)
      for(c=0;c<=2;c++){
	if(epsilon(a,b,c) == 0) continue;
	eabc = epsilon(a,b,c);
	for(ap=0;ap<=2;ap++){
	  for(bp=0;bp<=2;bp++)
	    for(cp=0;cp<=2;cp++){
	      if(epsilon(ap,bp,cp) == 0) continue;
	      eabcp = epsilon(ap,bp,cp);
	      corr1 += eabc * eabcp * inv2 * ((UDD1[alpha][a][beta][b][phix][cp]-UDD1[beta][a][alpha][b][phix][cp]) * (UDD2[delta][ap][epsil][bp][gammax][c] - UDD2[epsil][ap][delta][bp][gammax][c]));
	      corr2 += eabc * eabcp * inv2 * ((UDD1[delta][a][epsil][b][gammax][cp] - UDD1[epsil][a][delta][b][gammax][cp]) * (UDD2[alpha][ap][beta][bp][phix][c]-UDD2[beta][ap][alpha][bp][phix][c]));
	    };
	};
      };
  
  return corr1 - corr2;

}

