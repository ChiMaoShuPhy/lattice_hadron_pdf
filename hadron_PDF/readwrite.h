// -*- C++ -*-
// $Id$
/*! \file
 *  \brief Personal read/write routine definitions
 */

#include "chromabase.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <qdp.h>
using namespace QDP;

#ifndef __readwrite_h__
#define __readwrite_h__

class readwrite
{
 public:

  int writeCorrelator(const multi1d<DComplex>& corr, std::string name, int type, int Nt,  multi1d<int>& srce_pt);
  int writeCorrelator(const multi1d<DComplex>& corr, std::string name, int type, int Nt,  
                      multi1d<int>& srce_pt1, multi1d<int>& srce_pt2);
  int writeAvgCorrelator(const multi1d<DComplex>& corr, std::string name, int type, int Nt, int t0);

 private:

};


#endif
