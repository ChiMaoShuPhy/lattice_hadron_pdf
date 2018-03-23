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
#include "readwrite.h"

using namespace std;
using namespace QDP;

int readwrite::writeCorrelator(const multi1d<DComplex>& corr, std::string name, int type, int Nt,  multi1d<int>& srce_pt)
{
  fstream outre;
  fstream gplout;
  std::string typere,gpltypere;
  int tr;
  Real re,im;
 
  if(Layout::primaryNode()) {
    typere=name;
    gpltypere=name+"_gpl";
    outre.open(typere.c_str(),std::ios::out|std::ios::app);
    gplout.open(gpltypere.c_str(),std::ios::out|std::ios::app);
    outre << srce_pt[0];
    outre << " " << srce_pt[1];
    outre << " " << srce_pt[2];
    outre << " " << srce_pt[3];
    gplout << "# " << srce_pt[0] << " " << srce_pt[1] << " " << srce_pt[2] << " " << srce_pt[3] << "\n";
    for(int t=0;t<Nt;t++){
      if(t+srce_pt[3] >= Nt){
	tr=(t+srce_pt[3]-Nt);
	re = type * real(corr[tr]);
	im = type * imag(corr[tr]);
	outre << " " << re;
	gplout << t << " " << re << " " << im << "\n";
      }else{
	tr=t+srce_pt[3];
	re = real(corr[tr]);
	im = imag(corr[tr]);
	outre << " " << re;
	gplout << t << " " << re << " " << im << "\n";
      };
    };
    outre << "\n";
    outre << "       ";
    for(int t=0;t<Nt;t++){
      if(t+srce_pt[3] >= Nt){
	tr=(t+srce_pt[3]-Nt);
	im = type * imag(corr[tr]);
	outre << " " << im;
      }else{
	tr=t+srce_pt[3];
	im = imag(corr[tr]);
	outre << " " << im;
      };
    };
    outre << "\n";
    outre.close();

    gplout << "\n\n\n";
    gplout.close();
  }

  return 0;

}

int readwrite::writeCorrelator(const multi1d<DComplex>& corr, std::string name, int type, int Nt,  
                               multi1d<int>& srce_pt1, multi1d<int>& srce_pt2)
{
  fstream outre;
  fstream gplout;
  std::string typere,gpltypere;
  int tr;
  Real re,im;
 
  if(Layout::primaryNode()) {
    typere=name;
    gpltypere=name+"_gpl";
    outre.open(typere.c_str(),std::ios::out|std::ios::app);
    gplout.open(gpltypere.c_str(),std::ios::out|std::ios::app);
    outre << srce_pt1[0];
    outre << " " << srce_pt1[1];
    outre << " " << srce_pt1[2];
    outre << " " << srce_pt2[0];
    outre << " " << srce_pt2[1];
    outre << " " << srce_pt2[2];
    outre << " " << srce_pt2[3];
    gplout << "# " << srce_pt1[0] << " " << srce_pt1[1] << " " << srce_pt1[2] << " " << 
              srce_pt2[0] << " " << srce_pt2[1] << " " << srce_pt2[2] << " " << srce_pt2[3] << "\n";
    for(int t=0;t<Nt;t++){
      if(t+srce_pt1[3] >= Nt){
	tr=(t+srce_pt1[3]-Nt);
	re = type * real(corr[tr]);
	im = type * imag(corr[tr]);
	outre << " " << re;
	gplout << t << " " << re << " " << im << "\n";
      }else{
	tr=t+srce_pt1[3];
	re = real(corr[tr]);
	im = imag(corr[tr]);
	outre << " " << re;
	gplout << t << " " << re << " " << im << "\n";
      };
    };
    outre << "\n";
    outre << "             ";
    for(int t=0;t<Nt;t++){
      if(t+srce_pt1[3] >= Nt){
	tr=(t+srce_pt1[3]-Nt);
	im = type * imag(corr[tr]);
	outre << " " << im;
      }else{
	tr=t+srce_pt1[3];
	im = imag(corr[tr]);
	outre << " " << im;
      };
    };
    outre << "\n";
    outre.close();

    gplout << "\n\n\n";
    gplout.close();
  }

  return 0;

}

int readwrite::writeAvgCorrelator(const multi1d<DComplex>& corr, std::string name, int type, int Nt, int t0)
{
  fstream inAvg,outAvg;
  std::string typeAvg;
  int tr;
  int t;
  multi1d<DComplex> corrAvg(2*Nt);
  double avg;
  Real re;
  int Navg;

  if(Layout::primaryNode()) {
    typeAvg=name+"_avg";
    inAvg.open(typeAvg.c_str(),std::ios::in);

    if(! inAvg ) { /* file does not exist */
      inAvg.close();
      outAvg.open(typeAvg.c_str(),std::ios::out);
      
      outAvg << "1";
      
      for(t=0;t<Nt;t++){
	if(t+t0 >= Nt){
	  tr=(t+t0-Nt);
	  re = type * real(corr[tr]);
	  outAvg << " " << re;
	}else{
	  tr=t+t0;
	  re = real(corr[tr]);
	  outAvg << " " << re;
	};
      };
      outAvg << "\n";

      outAvg << " ";
      for(t=0;t<Nt;t++){
	if(t+t0 >= Nt){
	  tr=(t+t0-Nt);
	  re = type * imag(corr[tr]);
	  outAvg << " " << re;
	}else{
	  tr=t+t0;
	  re = imag(corr[tr]);
	  outAvg << " " << re;
	};
      };
      outAvg << "\n";
      outAvg.close();
    } else { /* file exists ! */
      
      inAvg >> Navg;
      for(t=0;t<2*Nt;t++) {  /* Read in averaged correlator */
	inAvg >> avg;
	corrAvg[t] = avg;
	corrAvg[t] *= Navg;
      };
      
      inAvg.close();
       
      for(t=0;t<Nt;t++){ // Add new correlator data and re-average
	if(t+t0 >= Nt){  
	  tr=(t+t0-Nt);
	  corrAvg[t] +=  type * real(corr[tr]);
	  corrAvg[t] /= (Navg+1);
	  corrAvg[t+Nt] += type * imag(corr[tr]);
	  corrAvg[t+Nt] /= (Navg+1);
	}else{
	  tr=t+t0;
	  corrAvg[t] += real(corr[tr]);
	  corrAvg[t] /= (Navg+1);
	  corrAvg[t+Nt] += imag(corr[tr]);
	  corrAvg[t+Nt] /= (Navg+1);
	};
      };
      Navg += 1;
      
      outAvg.open(typeAvg.c_str());

      outAvg << Navg;
      
      for(t=0;t<Nt;t++) {// Now write out new averaged correlator
	re = real(corrAvg[t]);
	outAvg << " " << re;
      }
         
      outAvg << "\n ";
      for(t=Nt;t<2*Nt;t++) {// Now write out new averaged correlator
	re = real(corrAvg[t]);
	outAvg << " " << re;
      }
      outAvg << "\n";
      outAvg.close();
      
    };
  }

  return 0;
}

