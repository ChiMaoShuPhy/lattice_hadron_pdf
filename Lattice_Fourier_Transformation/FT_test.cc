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

LatticePropagator getPointSource(const multi1d<int> &pos, const int &time_dir)
{
    START_CODE();

    //result vector:
    LatticePropagator result = zero;
    for (unsigned int c = 0; c < Nc; c++)
    {
        //fill lattice color vector with a one at pos:
        LatticeColorVector tmpcolorvec = zero;
        srcfil(tmpcolorvec, pos, c);

        //inject in larger source:
        for (unsigned int s = 0; s < Ns; s++)
        {
            LatticeFermion tmpferm = zero;
            CvToFerm(tmpcolorvec, tmpferm, s);
            FermToProp(tmpferm, result, c, s);
        }
    }
    END_CODE();

    return result;
}

bool linkageHack(void)
{
    bool foo = true;
    // Inline Measurements
    foo &= InlineAggregateEnv::registerAll();
    foo &= GaugeInitEnv::registerAll();

    return foo;
}

/*struct TimeSliceFunc : public SetFunc
{
  TimeSliceFunc(int dir): mu(dir) {}
  // Simply return the mu’th coordinate
  int operator()(const multi1d<int>& coord)
  {return coord[mu];}
  // The number of subsets is the length of the lattice
  // in direction mu
  int numSubsets() {return Layout::lattSize()[mu];}
  int mu; // state
};
*/

int main(int argc, char **argv)
{

    initialize(&argc, &argv);
    START_CODE();

    linkageHack();

    Real one = 1.0;
    Real zero = 0.0;
    Complex cmplx_zero = cmplx(Real(0.0),Real(0.0));

    XMLReader xml_in;

    try
    {
        xml_in.open(getXMLInputFileName());
    }
    catch (...)
    {
        cerr << "Error in Reading input XML:" << endl;
        QDP_abort(1);
    }

    multi1d<int> nrow;
    multi1d<int> site;
    Cfg_t cfg;
    Real mass = 0.0;
    string prpgtr_fl;
    multi1d<int> bndry_cndtns;
    int ncg_had = 0;
    int j_decay = 0;
    multi1d<int> t_srce;

    int dir;
    int steps;

    Real wvf_param;
    int wvfIntPar;

    //================= set parameters =====================

    try
    {
        read(xml_in, "/chroma/Param/nrow", nrow);
        read(xml_in, "/chroma/Param/site", site);
        //    cout<<"read nrow done"<<endl;
        read(xml_in, "/chroma/Cfg", cfg);
        //    cout<<"read cfg done"<<endl;
        read(xml_in, "/chroma/Param/InlineMeasurements/elem/propagator/Param/FermionAction/Mass", mass);
        //    cout<<"read mass done"<<endl;
        read(xml_in, "/chroma/Param/InlineMeasurements/elem/propagator/Param/FermionAction/FermState/FermionBC/boundary", bndry_cndtns);
        //    cout<<"read boundry condition done"<<endl;
        read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Source_Smearing/j_decay", j_decay);
        read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Source_Smearing/t_srce", t_srce);
        //  W[z+Delta_z,z], Delta_z: steps
        read(xml_in, "/chroma/Param/InlineMeasurements/elem/displacement/dir", dir);
        read(xml_in, "/chroma/Param/InlineMeasurements/elem/displacement/steps", steps);
    }
    catch (const string &e)
    {
        QDPIO::cerr << "Parsing XML: " << e << endl;
        QDP_abort(1);
    }

    Layout::setLattSize(nrow);
    Layout::create();

    LatticePropagator quark_prpgtr = zero;
    LatticePropagator quark_src = zero;
    quark_src = getPointSource(t_srce, j_decay);

    XMLReader rd(xml_in, "/chroma/Param/InlineMeasurements/elem/propagator/Param");

    SimpleFermBCParams bc_prmtrs;

    GroupXML_t InvertParam = readXMLGroup(rd, "InvertParam", "invType");

   cout<<"====="<<Nd<<"====="<<endl;

    //boundary conditions
    bc_prmtrs.boundary.resize(Nd);
    for (int i = 0; i < Nd - 1; i++)
    {
        bc_prmtrs.boundary[i] = bndry_cndtns[i]; //BC_TYPE_PERIODIC;
    }
    bc_prmtrs.boundary[Nd - 1] = bndry_cndtns[Nd - 1]; //BC_TYPE_ANTIPERIODIC;

    //calculate propagator

    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    Handle<FermBC<T, P, P>> fbc_handle(new SimpleFermBC<T, P, P>(bc_prmtrs));
    Handle<CreateFermState<T, P, P>> cfs(new CreateSimpleFermState<T, P, P>(fbc_handle));

    EvenOddPrecWilsonFermAct S(cfs, mass);

    //============= create gauge configurations ===================

    XMLReader gauge_file_xml, gauge_xml;
    multi1d<LatticeColorMatrix> u(Nd);
    //HotSt(u);

    gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);
    unitarityCheck(u);

    //======================================================
    //
    //  BEFORE rgauge
    //
    //======================================================

    //=================== 1st inversion =====================

    Handle<FermState<T, P, P>> state(S.createState(u));
    LinOpSysSolverCGEnv::registerAll();
    QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

    try
    {
        XMLFileWriter xml_out("1st_propagator");
        push(xml_out, "pion_crrltr");

        S.quarkProp(quark_prpgtr, xml_out,
                    quark_src, state, InvertParam,
                    quarkSpinType, ncg_had);

        xml_out.flush();
        pop(xml_out);
    }
    catch (const string &e)
    {
        QDPIO::cout << "Error trying to perform the inversion:" << e << endl;
        QDP_abort(1);
    }

    //calculate propagator with FH source
    LatticePropagator quark_prpgtr_Delta = zero;
    LatticePropagator FH_src = zero;
    LatticePropagator FH_quark_prpgtr = zero;
    LatticeColorMatrix gaugelink = zero;

    try
    {
        gaugelink = gauge_link(u, +1, dir, steps);

        quark_prpgtr_Delta = multi_shift<LatticePropagator>(quark_prpgtr, steps, dir, nrow[dir], bndry_cndtns[dir]);

        //Gamma(8)=gamma^3 comes from the definition of quasi-PDF

        //  LatticePropagator FH_src = Gamma(8)* gaugelink * quark_prpgtr_Delta; Wrong!, QDP does not support Gamma multiply with ColorMatrix

        FH_src = gaugelink * (Gamma(8) * quark_prpgtr_Delta);

        //======================Feynman-Hellman Source=========================

        Handle<FermState<T, P, P>> state(S.createState(u));

        LinOpSysSolverCGEnv::registerAll();
        
        QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

        //===============2nd inversion Feynman-Hellman propagator==================

        XMLFileWriter xml_out("fh_xmldat");
        push(xml_out, "FH_pion_crrltr");
        S.quarkProp(FH_quark_prpgtr, xml_out,
                    FH_src, state, InvertParam,
                    quarkSpinType, ncg_had);
        xml_out.flush();
        pop(xml_out);
    }
    catch (const string &e)
    {
        QDPIO::cout << "Error trying to perform the inversion:" << e << endl;
        QDP_abort(1);
    }

    //=========================pion-correlation function========================

    LatticePropagator anti_quark_prpgtr = Gamma(15) * quark_prpgtr * Gamma(15);
    LatticeComplex tr_prpgtrs = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * FH_quark_prpgtr);

    //=========================HDF5 Output========================

    LatticeDouble tr_prpgtrs_Re = real(tr_prpgtrs);
    LatticeDouble tr_prpgtrs_Im = imag(tr_prpgtrs);

    HDF5Writer fftw_test_file;
    fftw_test_file.open("fftw_test.h5", HDF5Base::ate);
    fftw_test_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.

    fftw_test_file.push("tr_prpgtrs");
    fftw_test_file.write("Re", tr_prpgtrs_Re, HDF5Base::trunc);
    fftw_test_file.write("Im", tr_prpgtrs_Im, HDF5Base::trunc);
    fftw_test_file.pop();

    fftw_test_file.close();
    
//    double test = peekSite(tr_prpgtrs_Re, site);

//    cout << "!!!!" << site[0] << " " << site[1]
//         << " " << site[2] << " " << site[3] << "!!!!" << endl;
//   cout << "!!!!" << peekSite(tr_prpgtrs_Re, site) << "!!!!" << endl;

//    cout << "!!!!" << typeid(peekSite(tr_prpgtrs_Re, site)).name() << "!!!!" << endl;
   //*cmplx(Real(1.0),Real(0.0))
    LatticeDouble r_dot_p = zero;
    LatticeDouble coord_in_dir;
    
    
    int length = nrow[j_decay]; 

   
    multi1d<Real>tr_prpgtrs_FT_Sum_Header_Re(length); 
    multi1d<Real>tr_prpgtrs_FT_Sum_Header_Im(length);
    multi1d<Real>tr_prpgtrs_Thorstn_fft_Re(length); 
    multi1d<Real>tr_prpgtrs_Thorstn_fft_Im(length);  

    multi1d<int>momentum2(4);
    momentum2[0] = 1; momentum2[1] = 2;
    momentum2[2] = 3; momentum2[3] = 4;


    cout<<"Mom = "<<momentum2[0]<<"  "<<momentum2[1]<<"  "<<momentum2[2]<<endl;
    
    
    StopWatch swatch;
    swatch.start();
    
    multi1d<DComplex> tr_prpgtrs_FT_Sum_Header(length); 
      tr_prpgtrs_FT_Sum_Header = timesilce_FT<multi1d<DComplex>,LatticeComplex>
                                        (tr_prpgtrs, nrow, momentum2, j_decay, Nd);

    swatch.stop();
    cout<<"my FT time= " << swatch.getTimeInSeconds() << " secs"<<endl;
   
    StopWatch swatchp;
    swatchp.start();

    Fourier Thorstn_fft(j_decay);
    LatticeComplex tr_prpgtrs_FTed_p = tr_prpgtrs;
    LatticeComplex tr_prpgtrs_FTed = tr_prpgtrs;

        multi1d<DComplex> Thorstn_ffted(length);


    LatticeDouble tr_prpgtrs_FTed_Re;
    LatticeDouble tr_prpgtrs_FTed_Im;

    tr_prpgtrs_FTed = Thorstn_fft(tr_prpgtrs_FTed_p,1);

    tr_prpgtrs_FTed_Re = real(tr_prpgtrs_FTed);
    tr_prpgtrs_FTed_Re = imag(tr_prpgtrs_FTed);


    multi1d<int>momentum3(4);
    momentum3[0] = 1; momentum3[1] = 2;
    momentum3[2] = 3; momentum3[3] = 4;

    for(int i = 0; i<length; i++)
    {
        momentum3[j_decay] = i;
        Thorstn_ffted[i] = peekSite(tr_prpgtrs_FTed,momentum3);
    }

    swatchp.stop();
    cout<<"my FT time= " << swatchp.getTimeInSeconds() << " secs"<<endl;


    for(int i = 0; i < length; i++)
    {
     tr_prpgtrs_FT_Sum_Header_Re[i] = real(tr_prpgtrs_FT_Sum_Header[i]);
     tr_prpgtrs_FT_Sum_Header_Im[i] = imag(tr_prpgtrs_FT_Sum_Header[i]);
     tr_prpgtrs_Thorstn_fft_Re[i] = real(Thorstn_ffted[i]);
     tr_prpgtrs_Thorstn_fft_Im[i] = imag(Thorstn_ffted[i]);
    }

    HDF5Writer FT_factor_file;
    FT_factor_file.open("Fourier.h5", HDF5Base::ate);
    FT_factor_file.set_stripesize(1048576); 

    FT_factor_file.push("tr_prpgtrs_FT_Sum_Header");
    FT_factor_file.write("Re", tr_prpgtrs_FT_Sum_Header_Re, HDF5Base::trunc);
    FT_factor_file.write("Im", tr_prpgtrs_FT_Sum_Header_Im, HDF5Base::trunc);
    FT_factor_file.pop();

    FT_factor_file.push("tr_prpgtrs_Thorstn_fft");
    FT_factor_file.write("Re", tr_prpgtrs_Thorstn_fft_Re, HDF5Base::trunc);
    FT_factor_file.write("Im", tr_prpgtrs_Thorstn_fft_Im, HDF5Base::trunc);
    FT_factor_file.pop();

    FT_factor_file.push("tr_prpgtrs_Thorstn_fft_all");
    FT_factor_file.write("Re", tr_prpgtrs_FTed_Re, HDF5Base::trunc);
    FT_factor_file.write("Im", tr_prpgtrs_FTed_Im, HDF5Base::trunc);
    FT_factor_file.pop();

    FT_factor_file.close();

    END_CODE();
    finalize();
    exit(0);
}
