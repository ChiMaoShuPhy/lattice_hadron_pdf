// -*- C++ -*-
// $Id$
/*! \file
 *  \brief Kaon_Nucleon_system- definitions
 */

#include "chromabase.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <qdp.h>
#include "readwrite.h"
#include "meson.h"
#include "nucleon.h"

using namespace QDP;

#ifndef __kaon_nucleon_system_h__
#define __kaon_nucleon_system_h__

class Kaon_Nucleon_system
{
public:
  
  // Default constructor -- forward declaration
  Kaon_Nucleon_system() {
//    srce_ptK.resize(4);
//    srce_ptN.resize(4);
//    lat_size.resize(4);
  };

  // constructor with passed with all up and strange propagators and source positions
/*  Kaon_Nucleon_system(const LatticePropagator& propU_Nucleon,const LatticePropagator& propU_Kaon, const LatticePropagator& propS_Kaon, 
		      const multi1d<int>& srcN, const multi1d<int>& srcK,const multi1d<int>& latsize) {
    srce_ptK.resize(4);
    srce_ptN.resize(4);
    lat_size.resize(4);
 
    initialize(propU_Nucleon, propU_Kaon, propS_Kaon, srcN, srcK, latsize);
  };
*/
  // Kaon and Nucleon objects
  Meson_system kaon;
  Nucleon_system nucleon;

  // initialization routine
  int initialize(const multi1d<int>& srcN, const multi1d<int>& momN,
                 const multi1d<int>& srcK, const multi1d<int>& momK,
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
  multi1d<int> srce_ptK;
  multi1d<int> srce_ptN;
  multi1d<int> pK;
  multi1d<int> pN;

};


#endif
