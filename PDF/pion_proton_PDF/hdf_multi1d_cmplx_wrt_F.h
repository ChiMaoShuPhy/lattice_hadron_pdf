#ifndef HDF_MULTI1D_CMPLX_WRT_F_H
#define HDF_MULTI1D_CMPLX_WRT_F_H

#include <string>
#include "chroma.h"

using namespace std;

void hdf_multi1d_cmplx_wrt_F(HDF5Writer& wrtr,int lngth,multi1d<DComplex> dat)
{
      multi1d<Double> dat_Re(lngth);
      multi1d<Double> dat_Im(lngth);
    for (int i = 0; i < lngth; i++)
    {
      dat_Re[i] = real(dat[i]);
      dat_Im[i] = imag(dat[i]);
    }
    wrtr.write("Re", dat_Re, HDF5Base::trunc);
    wrtr.write("Im", dat_Im, HDF5Base::trunc);
}


void hdf_multi1d_cmplx_wrt_F(HDF5Writer& wrtr,string tag,int lngth,multi1d<DComplex> dat)
{
      multi1d<Double> dat_Re(lngth);
      multi1d<Double> dat_Im(lngth);
    for (int i = 0; i < lngth; i++)
    {
      dat_Re[i] = real(dat[i]);
      dat_Im[i] = imag(dat[i]);
    }
    wrtr.push(tag);
    wrtr.write("Re", dat_Re, HDF5Base::trunc);
    wrtr.write("Im", dat_Im, HDF5Base::trunc);
    wrtr.pop();
}

#endif
