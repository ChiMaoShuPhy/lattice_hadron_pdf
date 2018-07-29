#ifndef MMNTM_PRJCTN_H
#define MMNTM_PRJCTN_H

#include "fourier_cpu.h"
#include "chroma.h"

using namespace Chroma;

class FFT
{
  public:
    FFT(LatticeComplex data, int j_decay, int t_lngth)
    {
        L = t_lngth;
        t = j_decay;
        Fourier fft(t);
        data_fft_all = fft(data, 1);
    }
    ~FFT(){}

     multi1d<DComplex> prjctn(multi1d<int>& hdrn_mmtm)
    {
        multi1d<DComplex> data_fft(L);
        hdrn_mmtm_1 = hdrn_mmtm;
        for (int i = 0; i < L; i++)
        { 
            hdrn_mmtm_1[3] = i;
            data_fft[i] = peekSite(data_fft_all, hdrn_mmtm_1);
        }
        return data_fft;
    }

  private:
    int L;
    int t;
    multi1d<int> hdrn_mmtm_1;
    LatticeComplex data_fft_all;
};

#endif // !MMNTM_PRJCTN