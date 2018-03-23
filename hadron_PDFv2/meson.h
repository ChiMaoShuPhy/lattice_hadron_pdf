// -*- C++ -*-
// $Id$
/*! \file
 *  \brief Meson- definitions
 */

#include "chromabase.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <qdp.h>
#include "mom_project.h"
#include "readwrite.h"

using namespace QDP;

#ifndef __meson_system_h__
#define __meson_system_h__

class Meson_system
{
public:
  
  // Default constructor -- forward declaration
  Meson_system() {
  };  

  // US blocks
  multi5d<DComplex> block;
  multi1d<Propagator> Mblock;

  // initialization routine
  int initialize(const multi1d<int>& src,const multi1d<int>& mom, const multi1d<int>& latsize,
                               std::string directory, std::string base);
  
  int setSmatrices(const LatticePropagator& propQ);
  int setSmatrices(const LatticePropagator& propQ, const LatticePropagator& propQbar);
  int calcBlock(const LatticePropagator& propQ);
  int calcBlock(const LatticePropagator& propQ, const LatticePropagator& propQbar);
  
  int calcDirect(const LatticePropagator& propQ);
  int calcDirect(const LatticePropagator& propQ, const LatticePropagator& propQbar);
  int calcExchange(const LatticePropagator& propQ, const LatticePropagator& propQbar);

  // meson mass routine
  int calcMassFromBlocks();
  int writeMass();

  // consistency test routines 
  bool doChecks;
  
  // free memory
  int free();

  // meson correlator
  multi1d<DComplex> direct;
  
private:
    
  int quarkPerm;
    
  std::string dirname;
  std::string basename;

  multi1d<int> lat_size;

  // Lattice color matrices
  multi2d<LatticeColorMatrix> Su;
  multi2d<LatticeColorMatrix> Ss;
  
  // private routine to calculate meson block
  int calcMeson_blocks();

  // zero momentum container
  mom_project Pxyz;

  // read/write container
  readwrite rw;

  // where the meson source is located
  multi1d<int> srce_pt;
  multi1d<int> p_mom;

};


#endif
