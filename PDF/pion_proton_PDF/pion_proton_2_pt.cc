//two-point function only, search for value of smear parameters
#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "barspinmat_w.cc"
#include "meas/smear/quark_smearing.h"
#include "meas/smear/gaus_smear.h"
#include "fourier_cpu.h"
#include "hdf_multi1d_cmplx_wrt_F.h"
#include "mmntm_prjctn.h"

#include "sobol.h"
#include "sobol.cpp"
#include <vector>

using namespace std;
using namespace QDP;
using namespace Chroma;
using namespace BaryonSpinMats;

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl
#define PASS_LN cout << "====== PASS  " << __LINE__ << " ======" << endl

typedef OLattice<PSpinMatrix<PColorMatrix<RScalar<double>, Nc>, Ns>> LatticePropagatorReImPrt;
typedef OLattice<PScalar<PColorMatrix<RScalar<double>, Nc>>> LatticeColorMatrixReImPrt;

struct int4
{
  int int_arry[4];
};

struct Real4
{
  DOUBLE double_arry[4];
};

//point source
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
// sobol sequence as source position
void sobol_src_gnrtr(int srnd_seed, int t_lngth, multi1d<int> nrow,
                     vector<int4> &src_vct)
{
  srand(srnd_seed);
  int seed = rand() % 256;
  int skip = 0;
  int4 src_tmp;

  src_tmp.int_arry[0] = 0;
  src_tmp.int_arry[1] = 0;
  src_tmp.int_arry[2] = 0;
  src_tmp.int_arry[3] = 0;

  for (unsigned long long i = 0; i < t_lngth + seed; ++i)
  {
    // Print a few dimensions of each point.
 //!!! bounded for 4-dimension
    for (unsigned d = 0; d < 4; ++d)
    {
      const int s = (int)(nrow[d] * sobol::sample(i, d));
      if (skip < seed)
      {        ;      }
      else
      {
        src_tmp.int_arry[d] = s;
      }
    }
    if (skip < seed)
    { ; }
    else
    {
      src_vct.push_back(src_tmp);
    }
    skip++;
  }
}

//////////////////////MAIN///////////////////////

int main(int argc, char **argv)
{

  initialize(&argc, &argv);
  START_CODE();

  linkageHack();

  Real one = 1.0;
  Real zero = 0.0;

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
  Cfg_t cfg;
  string cfg_file_name;
  Real mass = 0.0;
  string prpgtr_fl;
  multi1d<int> bndry_cndtns;
  int ncg_had = 0;
  int j_decay = 0;
  int Sobol_src_seed = 0;  //seed for sobol sequence generator
  int Sobol_src_Lngth = 0; //Length of 4-D sobol sequence as source position

  Real wvf_param;
  int wvfIntPar;

  //================= read in  parameters ===============
  try
  {
    read(xml_in, "/chroma/Param/nrow", nrow);
    //    cout<<"read nrow done"<<endl;
    read(xml_in, "/chroma/Cfg", cfg);
    read(xml_in, "/chroma/Cfg/cfg_file", cfg_file_name);
    cfg_file_name = cfg_file_name.substr(cfg_file_name.rfind("/") + 1, cfg_file_name.length());
    //    cout<<"read cfg done"<<endl;
    read(xml_in, "/chroma/Param/InlineMeasurements/elem/propagator/Param/FermionAction/Mass", mass);
    //    cout<<"read mass done"<<endl;
    read(xml_in, "/chroma/Param/InlineMeasurements/elem/propagator/Param/FermionAction/FermState/FermionBC/boundary", bndry_cndtns);
    //    cout<<"read boundry condition done"<<endl;
    read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Source_Smearing/j_decay", j_decay);
    read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Source_Smearing/Sobol_src_seed", Sobol_src_seed);
    read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Source_Smearing/Sobol_src_Lngth", Sobol_src_Lngth);
  }
  catch (const string &e)
  {
    QDPIO::cerr << "Parsing XML: " << e << endl;
    QDP_abort(1);
  }

  Layout::setLattSize(nrow);
  Layout::create();

  int t_lngth = nrow[j_decay];

  LatticePropagator quark_prpgtr = zero;
  LatticePropagator quark_src = zero;

  XMLReader rd(xml_in, "/chroma/Param/InlineMeasurements/elem/propagator/Param");

  SimpleFermBCParams bc_prmtrs;

  GroupXML_t InvertParam = readXMLGroup(rd, "InvertParam", "invType");

  //boundary conditions
  bc_prmtrs.boundary.resize(Nd);
  for (int i = 0; i < Nd - 1; i++)
  {
    bc_prmtrs.boundary[i] = bndry_cndtns[i]; //BC_TYPE_PERIODIC;
  }
  bc_prmtrs.boundary[Nd - 1] = bndry_cndtns[Nd - 1]; //BC_TYPE_ANTIPERIODIC;

  //create fermion state

  typedef LatticeFermion T;
  typedef multi1d<LatticeColorMatrix> P;
  Handle<FermBC<T, P, P>> fbc_handle(new SimpleFermBC<T, P, P>(bc_prmtrs));
  Handle<CreateFermState<T, P, P>> cfs(new CreateSimpleFermState<T, P, P>(fbc_handle));

  EvenOddPrecWilsonFermAct S(cfs, mass);

  //============= create gauge configurations =============

  XMLReader gauge_file_xml, gauge_xml;
  multi1d<LatticeColorMatrix> u(Nd);
  //HotSt(u);

  StopWatch swatch;

  swatch.start();

  gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);
  unitarityCheck(u);

  //=============  rgauge  ==============;
  //   rgauge(u);
  
  LinOpSysSolverCGEnv::registerAll();
  QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

  // generating sobol sequence as source position//

  multi1d<int>t_src(Nd);

  vector<int4> sobol_src;
  sobol_src_gnrtr(Sobol_src_seed, Sobol_src_Lngth, nrow, sobol_src);
 //loop over source position
  
    HDF5Writer pion_file;
    pion_file.open("pion_2pt.h5", HDF5Base::ate);
    pion_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.

    HDF5Writer proton_file;
    proton_file.open("proton_2pt.h5", HDF5Base::ate);
    proton_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
 //--1push
    pion_file.push(cfg_file_name);
    proton_file.push(cfg_file_name);


    XMLReader smr_xml_in;

    try
  {
    smr_xml_in.open("snk_smr.txt");
  }
  catch (...)
  {
    cerr << "Error in Reading smear XML:" << endl;
    QDP_abort(1);
  }
  
  int nmbr_snks = 0;
  REAL smr_wdth = 0.0;
  int ItrGaus = 0;

        for (int i = 0; i < Nd; i++)
    {
 //     t_src[i] = sobol_src[src_iter].int_arry[i];
        t_src[i] = 0;
    }
//    multi1d<Handle< QuarkSourceSink<LatticePropagator> > > Sink_Smr(smr_wdth.size())
    try
    {
     read(smr_xml_in, "/sinks/elem/nmbr_snks", nmbr_snks);
    }
    catch (...)
    {
      cerr << "Error in Reading in nmbr_snks:"<<endl;
      QDP_abort(1);
    }
//  for (int src_iter = 0; src_iter < sobol_src.size(); src_iter++)
   for (int smr_itr = 0; smr_itr < nmbr_snks; smr_itr++)
  {

 //============= 1st inversion ===============
  
  Handle<FermState<T, P, P>> state(S.createState(u));
    quark_src = getPointSource(t_src, j_decay);
  
    string smr_path="/sinks/elem/SinkSmearing_"+to_string(smr_itr);

     try
    {
     read(smr_xml_in, smr_path+"/SmearingParam/wvf_param",smr_wdth );
     read(smr_xml_in, smr_path+"/SmearingParam/wvfIntPar",ItrGaus );;
    }
    catch (...)
    {
      cerr << "Error in Reading in smearing parameters:"<<endl;
      QDP_abort(1);
    }
  

  //!!!!!!source smearing
    gausSmear(u, quark_src, smr_wdth, ItrGaus, j_decay);
    try
    {
      XMLFileWriter xml_out("1st_propagator");
      push(xml_out, "proton_crrltr");

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

  LatticePropagator quark_prpgtr_snk_smr = quark_prpgtr;

  //!!!!sink smearing
  QuarkSinkSmearingEnv::registerAll();
  
  cout<<smr_path<<endl;
    
  try
    {
      // Get the name of the smearing
      string smearing_name;
      read(smr_xml_in, "/sinks/elem/SinkSmearing_1/SinkType", smearing_name);

      QDPIO::cout << smearing_name<<endl;

      Handle< QuarkSourceSink<LatticePropagator> >
	    sinkSmearing(ThePropSinkSmearingFactory::Instance().createObject(smearing_name,
							 smr_xml_in, smr_path, u));

     (*sinkSmearing)(quark_prpgtr_snk_smr);
    }
  catch(const std::string& e)
    {
      QDPIO::cerr << ": Caught Exception creating sink: " << e << endl;
      QDP_abort(1);
    }

    string smr_prmtrs = to_string(smr_wdth) + " " + to_string(ItrGaus);

    cout<<smr_prmtrs<<endl;
 //push smear parameters
 //--2push    
    pion_file.push(smr_prmtrs);
    proton_file.push(smr_prmtrs);       
          

          swatch.stop();
          
          SpinMatrix cg5P_pstv = Cg5();

          //!!!! IMPORTNT
          SpinMatrix T_unpol = Tunpol();
          //!!!! IMPORTNT

          LatticePropagator anti_quark_prpgtr = Gamma(15) * quark_prpgtr * Gamma(15);
          LatticePropagator anti_quark_prpgtr_snk_smr = Gamma(15) * quark_prpgtr_snk_smr * Gamma(15);


            //pion-pion two point correlation
            LatticeComplex pion_TwoPnt_pnt_snk = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * quark_prpgtr);

            LatticeComplex pion_TwoPnt_snk_smr = trace(Gamma(15) * adj(anti_quark_prpgtr_snk_smr) * Gamma(15) * quark_prpgtr_snk_smr);

            //proton-prton two point correlation
            LatticeComplex proton_TwoPnt_pnt_snk_1 =
                trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

            LatticeComplex proton_TwoPnt_pnt_snk_2 =
                trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

            LatticeComplex proton_TwoPnt_snk_smr_1 =
                trace(T_unpol * traceColor(quark_prpgtr_snk_smr * traceSpin(quarkContract13(quark_prpgtr_snk_smr * cg5P_pstv, cg5P_pstv * quark_prpgtr_snk_smr))));

            LatticeComplex proton_TwoPnt_snk_smr_2 =
                trace(T_unpol * traceColor(quark_prpgtr_snk_smr * quarkContract13(quark_prpgtr_snk_smr * cg5P_pstv, cg5P_pstv * quark_prpgtr_snk_smr)));

            FFT pion_TwoPnt_pnt_snk_fft(pion_TwoPnt_pnt_snk, j_decay, t_lngth);
            FFT proton_TwoPnt_pnt_snk_1_fft(proton_TwoPnt_pnt_snk_1, j_decay, t_lngth);
            FFT proton_TwoPnt_pnt_snk_2_fft(proton_TwoPnt_pnt_snk_2, j_decay, t_lngth);

            FFT pion_TwoPnt_snk_smr_fft(pion_TwoPnt_snk_smr, j_decay, t_lngth);
            FFT proton_TwoPnt_snk_smr_1_fft(proton_TwoPnt_snk_smr_1, j_decay, t_lngth);
            FFT proton_TwoPnt_snk_smr_2_fft(proton_TwoPnt_snk_smr_2, j_decay, t_lngth);

            multi1d<DComplex> pion_TwoPnt_pnt_snk_fftd;
            multi1d<DComplex> pion_TwoPnt_snk_smr_fftd;

            multi1d<DComplex> proton_TwoPnt_pnt_snk_1_fftd;
            multi1d<DComplex> proton_TwoPnt_pnt_snk_2_fftd; 

            multi1d<DComplex> proton_TwoPnt_snk_smr_1_fftd;
            multi1d<DComplex> proton_TwoPnt_snk_smr_2_fftd;          

            multi1d<int> hdrn_mmntm_zero(4);
            hdrn_mmntm_zero[0] = 0; hdrn_mmntm_zero[1] = 0; 
            hdrn_mmntm_zero[2] = 0; hdrn_mmntm_zero[3] = 0; 

            pion_TwoPnt_pnt_snk_fftd = pion_TwoPnt_pnt_snk_fft.prjctn(hdrn_mmntm_zero);
            pion_TwoPnt_snk_smr_fftd = pion_TwoPnt_snk_smr_fft.prjctn(hdrn_mmntm_zero);

            proton_TwoPnt_pnt_snk_1_fftd = proton_TwoPnt_pnt_snk_1_fft.prjctn(hdrn_mmntm_zero);
            proton_TwoPnt_pnt_snk_2_fftd = proton_TwoPnt_pnt_snk_2_fft.prjctn(hdrn_mmntm_zero);

            proton_TwoPnt_snk_smr_1_fftd = proton_TwoPnt_snk_smr_1_fft.prjctn(hdrn_mmntm_zero);
            proton_TwoPnt_snk_smr_2_fftd = proton_TwoPnt_snk_smr_2_fft.prjctn(hdrn_mmntm_zero);
    
            hdf_multi1d_cmplx_wrt_F(pion_file,"TwoPnt_pnt_snk",t_lngth,pion_TwoPnt_pnt_snk_fftd);
            hdf_multi1d_cmplx_wrt_F(pion_file,"TwoPnt_snk_smr",t_lngth,pion_TwoPnt_snk_smr_fftd);
             
//--3push  
            
            hdf_multi1d_cmplx_wrt_F(proton_file,"TwoPnt_pnt_snk_1",t_lngth,proton_TwoPnt_pnt_snk_1_fftd);
            
            hdf_multi1d_cmplx_wrt_F(proton_file,"TwoPnt_pnt_snk_2",t_lngth,proton_TwoPnt_pnt_snk_2_fftd);

            hdf_multi1d_cmplx_wrt_F(proton_file,"TwoPnt_snk_smr_1",t_lngth,proton_TwoPnt_snk_smr_1_fftd);
            
            hdf_multi1d_cmplx_wrt_F(proton_file,"TwoPnt_snk_smr_2",t_lngth,proton_TwoPnt_snk_smr_2_fftd);
//-3pop               

          QDPIO::cout << "\n"
                      << "==============="
                      << "\n"
                      << "Cfg file: " << cfg_file_name << "\n"
                      << "Source position: " << t_src[0] << ", " << t_src[1] << ", " << t_src[2] << ", " << t_src[3] << "\n"
                      << "Correlation function calculation done, time= " << swatch.getTimeInSeconds() << " secs"
                      << "\n"
                      << "===============" << endl;
    //Fourier Transform, overall factor 1/sqrt(L_1*L_2*L3) contained in Fourier              multi1d<DComplex> pion_tr_prpgtrs_fft_all(t_lngth);
          
  //           for (int mom = 1; mom < 2; mom++)          
   //--4pop        
   
          cout << "===============\n"

               << "nrow = " << nrow[0] << " " << nrow[1]
               << " " << nrow[2] << " " << nrow[3] << "\n"

               << "bndry_cndtn = " << bndry_cndtns[0] << " " << bndry_cndtns[1]
               << " " << bndry_cndtns[2] << " " << bndry_cndtns[3] << "\n"

               << "Cfg file: " << cfg_file_name << "\n"
               << "Source position: " << t_src[0] << ", " << t_src[1] << ", " << t_src[2] << ", " << t_src[3] << "\n"

               << "smearing witdth, iterations: " << smr_wdth << ", " << ItrGaus <<"\n"

               << "===============" << endl;
  //--2pop
    pion_file.pop();   //pop smear parameters
    proton_file.pop(); //pop smear parameters
  }//loop over different  smear parameters
  //--1pop
  pion_file.pop();   //pop cfg file name
  proton_file.pop(); //pop cfg file name

  pion_file.close();
  proton_file.close();
  END_CODE();
  finalize();
  exit(0);
}