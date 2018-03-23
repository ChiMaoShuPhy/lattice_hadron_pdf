// -*- C++ -*-
// $Id$
/*! \file
 *  \brief Two_Nucleon_system- definitions
 */

#include "chromabase.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <qdp.h>
#include "readwrite.h"
#include "nucleon.h"

using namespace QDP;

#ifndef __two_nucleon_system_h__
#define __two_nucleon_system_h__

class Two_Nucleon_system
{
public:
  
  // Default constructor -- forward declaration
  Two_Nucleon_system() {
  };

  // Nucleon objects
  Nucleon_system nucleon1;
  Nucleon_system nucleon2;

  // initialization routine
  int initialize(const multi1d<int>& srcN1, const multi1d<int>& momN1,
                 const multi1d<int>& srcN2, const multi1d<int>& momN2,
                 const multi1d<int>& latsize, std::string dirname, std::string basename);


  // consistency test routines 
  bool doChecks;
  
  int calcExchange(const LatticePropagator& nucleonQ,const LatticePropagator& kaonQ,
                                    const LatticePropagator& kaonS);

  int calcDirect(const LatticePropagator& nucleonQ,const LatticePropagator& kaonQ,
                                    const LatticePropagator& kaonS);

  
  multi1d<DComplex> direct;
  
  // nucleon and meson exchange block
  multi1d<DComplex> exch1;
  multi1d<DComplex> exch2;
  multi1d<DComplex> exch3;
  
  int writeCorrelators();
  
  int free();
 
private:

  std::string dirname;
  std::string basename;
  
  multi1d<int> lat_size;

  // read/write container
  readwrite rw;

  // where the kaon source is located
  // and the nucleon source
  multi1d<int> srce_ptN1;
  multi1d<int> srce_ptN2;
  multi1d<int> pN1;
  multi1d<int> pN2;

};


#endif
