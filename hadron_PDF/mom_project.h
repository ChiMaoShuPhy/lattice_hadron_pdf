// -*- C++ -*-
// $Id$
/*! \file
 *  \brief Momentum projection routines definitions
 */

#include "chromabase.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <qdp.h>
using namespace QDP;

#ifndef __mom_project_h__
#define __mom_project_h__

class mom_project
{
 public:

  mom_project() {
    p.resize(Nd);
    src.resize(Nd);
  };

  mom_project( const multi1d<int>& momentum, const multi1d<int>& source_pts ) {
    p.resize(Nd);
    src.resize(Nd);
    initialize( momentum, source_pts );
  };

  multi1d<DComplex> project_momentum( const LatticeComplex& S);
  multi1d<ColorMatrix> project_momentum( const LatticeColorMatrix& S );
  multi1d<Propagator> project_momentum( const LatticePropagator& S );
  /* Propagator project_momentum( const LatticePropagator& S ); */

  int initialize( const multi1d<int>& momentum, const multi1d<int>& source_pts );

 private:

  multi1d<Int> p;
  multi1d<LatticeInteger> src;
  LatticeComplex phase;

  int calcPhase(void);
  
};


#endif
