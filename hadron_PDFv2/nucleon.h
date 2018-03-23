// -*- C++ -*-
// $Id$
/*! \file
 *  \brief Nucleon- definitions
 */

#include "chromabase.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <qdp.h>
#include "mom_project.h"
#include "udd.h"
#include "readwrite.h"
//#include "spin_matrix.h"

using namespace QDP;

#ifndef __nucleon_system_h__
#define __nucleon_system_h__

class Nucleon_system
{
public:
  
  // Default constructor -- forward declaration
  Nucleon_system() {
  };
  
  // UDD blocks
  udd<DComplex> blockG1pU; //G1p_e1_r1;
  udd<DComplex> blockG1pD; //G1p_e1_r2;

  // initialization routine
  int initialize(const multi1d<int>& src,const multi1d<int>& mom, const multi1d<int>& latsize,std::string dirname, std::string basename, std::string spin);
  
  int setSmatrices(const LatticePropagator& prop);
  int setSmatrices(const LatticePropagator& prop1, const LatticePropagator& prop2);
  int setSmatrices(const LatticePropagator& prop1, const LatticePropagator& prop2, const LatticePropagator& prop3);
  
  int calcBlock(const LatticePropagator& prop);
  int calcBlock(const LatticePropagator& prop1, const LatticePropagator& prop2, int qperm);
  int calcBlock(const LatticePropagator& prop1, const LatticePropagator& prop2, const LatticePropagator& prop3);
  
  // routines for calculating terms needed for dressed quark propagators
  int calcDirect(const LatticePropagator& prop1);
  int calcExchange(const LatticePropagator& prop1, const LatticePropagator& prop2);

  // nucleon mass routines
  int calcMassFromBlocks();
  int writeMass();

  // calculates NN correlators
  int calc_G1pG1p_A1p();
  
  // free memory
  int free();
  
  bool doChecks;

  // nucleon correlator
  multi1d<DComplex> corrU;
  multi1d<DComplex> corrD;
  multi1d<DComplex> direct;  // spin averaged correlator
  
  multi1d<Propagator> exch1,exch2,exch3; // spin averaged exchange terms
  
//  LatticePropagator QPP; // dressed up quark with proton
//  LatticePropagator QNN; // dressed up quark with neutron
//  LatticePropagator QPN; // dressed exchange quark with proton/neutron
  
private:

  // where to write things out. . .
  std::string dirname;
  std::string basename;
    
  multi1d<int> lat_size;

  // Lattice color matrices
  multi2d<LatticeColorMatrix> S1;
  multi2d<LatticeColorMatrix> S2;
  multi2d<LatticeColorMatrix> S3;

  multi1d<DComplex> UDD_permuted(int alpha, int a, int beta, int b, int gamma, int c, int alphap, int betap, int gammap);
  LatticeComplex UDD_unpermuted(int alpha, int a, int beta, int b, int gamma, int c, int alphap, int betap, int gammap);
  multi1d<DComplex> UDD_UDD_permuted(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
				     const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);

    // calculates direct diagram
  multi1d<DComplex> calcDirect(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
			       const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);

  // calculates exchange diagram
  //  multi1d<DComplex> calcExchange(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
  //				 const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);

  // calculates exchange diagram
  multi1d<DComplex> calcExchangeDD1(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
				 const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);

  // calculates exchange diagram
  multi1d<DComplex> calcExchangeDD2(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
				 const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);

  // calculates exchange diagram
  multi1d<DComplex> calcExchangeDD3(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
				 const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);

  // calculates exchange diagram
  multi1d<DComplex> calcExchangeDD4(int alpha, int beta, int gamma, int delta, int epsilon, int phi, 
				 const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);


  // calculates exchange diagram
  multi1d<DComplex> calcExchangeU(int XX1, int beta, int gamma, int XX2, int epsilon, int phi, 
				 const udd<DComplex>& UDD1, const udd<DComplex>& UDD2);

  // zero momentum container
  mom_project Pxyz;

  // read/write container
  readwrite rw;

  // where the source is located
  multi1d<int> srce_pt;
  multi1d<int> p_mom;

  // calculates the UDD blocks
  int calcNucleon_blocks(udd<DComplex>& block, int mu1, int mu2, int mu3);

  // multiplies a scalar with a multi1d color matrix
  multi1d<DComplex> scalarMult(double coeff, multi1d<DComplex> md1cm);
  
    // tells how to order the quark propagators for different sources. . .
  int quarkPerm;  // 111 all quarks at the same source
                  // 112 first two quarks at the same source, the last different
                  // . . .
                  // 123 all quarks have different sources.
  
  LatticeSpinMatrixD ProjG;
  LatticeSpinMatrixD ProjS;

};


#endif
