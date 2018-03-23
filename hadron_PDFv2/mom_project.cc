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
#include "chroma.h"
#include "mom_project.h"
using namespace QDP;

int mom_project::initialize( const multi1d<int>& momentum, const multi1d<int>& source_pts )
{
  int i;

  for(i=0;i<Nd;i++) {
    p[i]=momentum[i];
    src[i]=source_pts[i];
  };

  calcPhase();

  return 0;
}

int mom_project::calcPhase()
{
  const Real twopi = 6.283185307179586476925286;
  multi1d<Int> mom(Nd);
  multi1d<LatticeInteger> my_coord(Nd);
  int mu;
  LatticeReal p_dot_x = zero;
  LatticeColorMatrix fwd;

  for(mu=0; mu<Nd; ++mu){
    my_coord[mu]=Layout::latticeCoordinate(mu);
  };

  for(mu=0; mu < Nd; ++mu){
    if (mu == Nd - 1) continue;
    p_dot_x += LatticeReal(my_coord[mu] - src[mu])*twopi*Real(p[mu]) / Layout::lattSize()[mu];
  };

  phase = cmplx(cos(p_dot_x),sin(p_dot_x));

  return 0;

}

multi1d<Propagator> mom_project::project_momentum( const LatticePropagator& S )
{
  multi1d<Propagator> CC;

  CC=sumMulti(phase * S,TimeSliceFunc(Nd-1));

  return CC;
}

multi1d<ColorMatrix> mom_project::project_momentum( const LatticeColorMatrix& S )
{
  multi1d<ColorMatrix> CC;

  CC=sumMulti(phase * S,TimeSliceFunc(Nd-1));

  return CC;
}

multi1d<DComplex> mom_project::project_momentum( const LatticeComplex& S )
{
  multi1d<DComplex> CC;

  CC=sumMulti(phase * S,TimeSliceFunc(Nd-1));

  return CC;
}
