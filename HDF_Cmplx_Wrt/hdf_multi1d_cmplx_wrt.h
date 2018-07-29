#ifndef HDF_MULTI1D_CMPLX_WRT_H
#define HDF_MULTI1D_CMPLX_WRT_H

#include <string>
#include "chroma.h"

class hdf_multi1d_cmplx_wrt
{
public:
    HDF5Writer writer;
    int L;

    hdf_multi1d_cmplx_wrt(HDF5Writer& wrtr, int lngth)
  {
    L= lngth;
    writer = wrtr;
  }
  
  void wrt(string tag, multi1d<DComplex>& dat)
  {
    multi1d<Double>dat_Re(L);
    multi1d<Double>dat_Im(L);
    for (int i = 0; i < L; i++)
    {
      dat_Re[i] = real(dat[i]);
      dat_Im[i] = imag(dat[i]);
    }
    writer.push(tag);
    writer.write("Re", dat_Re, HDF5Base::trunc);
    writer.write("Im", dat_Im, HDF5Base::trunc);
    writer.pop();
  }

    void wrt2(multi1d<DComplex>& dat)
  {
    multi1d<Double>dat_Re(L);
    multi1d<Double>dat_Im(L);
    for (int i = 0; i < L; i++)
    {
      dat_Re[i] = real(dat[i]);
      dat_Im[i] = imag(dat[i]);
    }
    writer.write("Re", dat_Re, HDF5Base::trunc);
    writer.write("Im", dat_Im, HDF5Base::trunc);
  }

//  multi1d<Double>dat_Re;
//  multi1d<Double>dat_Im;
};

#endif
