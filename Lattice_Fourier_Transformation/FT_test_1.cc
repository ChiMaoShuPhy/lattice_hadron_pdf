#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "gauge_link.h"
#include "multi_shift.h"
#include "timeslice_FT.h"
#include "fourier_cpu.h"
//#include <complex.h>
//#include <fftw3.h>
//#include <typeinfo>

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl

using namespace std;
using namespace QDP;
using namespace Chroma;

bool linkageHack(void)
{
    bool foo = true;
    // Inline Measurements
    foo &= InlineAggregateEnv::registerAll();
    foo &= GaugeInitEnv::registerAll();

    return foo;
}

int main(int argc, char **argv)
{

    initialize(&argc, &argv);
    START_CODE();

    linkageHack();
//lattice layout setup
    multi1d<int> nrow(4);

    nrow[0] = 4;nrow[1] = 4;nrow[2] = 4;nrow[3] = 32;
    int j_decay = 3;
    int length = nrow[j_decay];
    int Nd = 4;

    Layout::setLattSize(nrow);
    Layout::create();
//use lattice cooordinate in time direction as test data
    LatticeComplex tr_prpgtrs = (LatticeComplex)Layout::latticeCoordinate(j_decay);

    LatticeDouble tr_prpgtrs_Re = real(tr_prpgtrs);
    LatticeDouble tr_prpgtrs_Im = imag(tr_prpgtrs);
//my FT
    multi1d<Real> tr_prpgtrs_timeslice_FT_Re(length);
    multi1d<Real> tr_prpgtrs_timeslice_FT_Im(length);
//Thorsten's FT
    multi1d<Real> tr_prpgtrs_Thorsten_FT_Re(length);
    multi1d<Real> tr_prpgtrs_Thorsten_FT_Im(length);
//momentum 
    multi1d<int> momentum2(4);
    momentum2[0] = 0; momentum2[1] = 0; momentum2[2] = 2; momentum2[3] = 0;

    multi1d<int> momentum3(4);
    momentum3[0] = 0; momentum3[1] = 0; momentum3[2] = 2; momentum3[3] = 0;
//perfomr my FT
    multi1d<DComplex> tr_prpgtrs_timeslice_FT(length);
    tr_prpgtrs_timeslice_FT = timesilce_FT<multi1d<DComplex>, LatticeComplex>(tr_prpgtrs, nrow, momentum2, j_decay, Nd);
//perform Thorsten's FT
    Fourier Thorsten_fft(j_decay);
    LatticeComplex tr_prpgtrs_Thorsten_FT_All = tr_prpgtrs;
    multi1d<DComplex> tr_prpgtrs_Thorsten_FT(length);

        LatticeComplex tr_prpgtrs_Thorsten_FT_All_p = 
Thorsten_fft(tr_prpgtrs_Thorsten_FT_All, 1);

    LatticeDouble tr_prpgtrs_Thorsten_FT_All_Re = real(tr_prpgtrs_Thorsten_FT_All_p);
    LatticeDouble tr_prpgtrs_Thorsten_FT_All_Im = imag(tr_prpgtrs_Thorsten_FT_All_p);

//extarct the momentum we want and make time slice
    for (int i = 0; i < length; i++)
    {
        momentum3[j_decay] = i;
        tr_prpgtrs_Thorsten_FT[i] = peekSite(tr_prpgtrs_Thorsten_FT_All_p, momentum3);
    }

//take the Re/Im part of my/Thorsten's FT
    for (int i = 0; i < length; i++)
    {
        tr_prpgtrs_timeslice_FT_Re[i] = real(tr_prpgtrs_timeslice_FT[i]);
        tr_prpgtrs_timeslice_FT_Im[i] = imag(tr_prpgtrs_timeslice_FT[i]);
        tr_prpgtrs_Thorsten_FT_Re[i] = real(tr_prpgtrs_Thorsten_FT[i]);
        tr_prpgtrs_Thorsten_FT_Im[i] = imag(tr_prpgtrs_Thorsten_FT[i]);
    }
//print out momentum
    cout<<momentum2[0]<<", "<<momentum2[1]<<", "<<momentum2[2]<<endl;
    cout<<momentum3[0]<<", "<<momentum3[1]<<", "<<momentum3[2]<<endl;

    HDF5Writer FT_factor_file;
    FT_factor_file.open("Fourier_1.h5", HDF5Base::ate);
    FT_factor_file.set_stripesize(1048576);
//write out the test data
    FT_factor_file.push("tr_prpgtrs");
    FT_factor_file.write("Re", tr_prpgtrs_Re, HDF5Base::trunc);
    FT_factor_file.write("Im", tr_prpgtrs_Im, HDF5Base::trunc);
    FT_factor_file.pop();
//write out the my FT results
    FT_factor_file.push("tr_prpgtrs_timeslice_FT");
    FT_factor_file.write("Re", tr_prpgtrs_timeslice_FT_Re, HDF5Base::trunc);
    FT_factor_file.write("Im", tr_prpgtrs_timeslice_FT_Im, HDF5Base::trunc);
    FT_factor_file.pop();
//write out the Thorsten's FT results, momentum extracted
    FT_factor_file.push("tr_prpgtrs_Thorstn_fft");
    FT_factor_file.write("Re", tr_prpgtrs_Thorsten_FT_Re, HDF5Base::trunc);
    FT_factor_file.write("Im", tr_prpgtrs_Thorsten_FT_Im, HDF5Base::trunc);
    FT_factor_file.pop();
//write out the Thorsten's FT results, all momenta
    FT_factor_file.push("tr_prpgtrs_Thorstn_fft_All");
    FT_factor_file.write("Re", tr_prpgtrs_Thorsten_FT_All_Re, HDF5Base::trunc);
    FT_factor_file.write("Im", tr_prpgtrs_Thorsten_FT_All_Im, HDF5Base::trunc);
    FT_factor_file.pop();


    FT_factor_file.close();

    END_CODE();
    finalize();
    exit(0);
}
