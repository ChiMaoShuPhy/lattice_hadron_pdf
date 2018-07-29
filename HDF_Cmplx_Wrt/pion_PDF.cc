#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "gauge_link.h"
#include "multi_shift.h"
#include "hdf_multi1d_cmplx_wrt.h"
#include "hdf_multi1d_cmplx_wrt_F_bak.h"
#include "fourier_cpu.h"
#include "mmntm_prjctn.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl
#define PASS_ln cout << "====== PASS  " << __LINE__ << " ======" << endl


typedef OLattice<PSpinMatrix<PColorMatrix<RScalar<double>, Nc>, Ns>> LatticePropagatorReImPrt;

typedef OLattice<PScalar<PColorMatrix<RScalar<double>, Nc>>> LatticeColorMatrixReImPrt;

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

  LatticePropagator quark_prpgtr_Delta = zero;
  LatticePropagator FH_src = zero;
  LatticePropagator FH_quark_prpgtr = zero;
  LatticeColorMatrix gaugelink = zero;

  try
  {
    gaugelink = gauge_link(u, +1, dir, steps);

    quark_prpgtr_Delta = multi_shift<LatticePropagator>(quark_prpgtr, steps, dir, nrow[dir], bndry_cndtns[dir]);

    
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

  Fourier fft(j_decay);
  int lngth = nrow[j_decay];

  LatticeComplex tr_prpgtrs_fft_all = fft(tr_prpgtrs, 1);
  multi1d<int> hdrn_mmntm_zero(4);
  multi1d<DComplex> tr_prpgtrs_fft(lngth);
  multi1d<Double> tr_prpgtrs_fft_Re(lngth);
  multi1d<Double> tr_prpgtrs_fft_Im(lngth);

  hdrn_mmntm_zero[0] = 0;
  hdrn_mmntm_zero[1] = 0;
  hdrn_mmntm_zero[2] = 0;
  hdrn_mmntm_zero[3] = 0;

  for (int i = 0; i < lngth; i++)
  {
    hdrn_mmntm_zero[3] = i;

    tr_prpgtrs_fft[i] = peekSite(tr_prpgtrs_fft_all, hdrn_mmntm_zero);

    tr_prpgtrs_fft_Re[i] = real(tr_prpgtrs_fft[i]);
    tr_prpgtrs_fft_Im[i] = imag(tr_prpgtrs_fft[i]);
  }

  multi1d<DComplex> tr_prpgtrs_fft_class;

  FFT FFT_class_test(tr_prpgtrs,j_decay, lngth);
  tr_prpgtrs_fft_class = FFT_class_test.prjctn(hdrn_mmntm_zero);


  string pi_mtrx_elmnt = "<pi(y,t)|q-bar(x).ga_3.W[x,x+" + to_string(steps) + "_stps_in_" + to_string(dir) + "].q(x+" + to_string(steps) + "_stps_in_" + to_string(dir) + ")|pi(0,0)>";

  HDF5Writer pion_PDF_file;
  pion_PDF_file.open("pion_PDF.h5", HDF5Base::ate);
  pion_PDF_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
PASS_ln;
  pion_PDF_file.push("tr_prpgtrs_fft");
  pion_PDF_file.write("Re", tr_prpgtrs_fft_Re, HDF5Base::trunc);
  pion_PDF_file.write("Im", tr_prpgtrs_fft_Im, HDF5Base::trunc);
  pion_PDF_file.pop();
PASS_ln;
  pion_PDF_file.push("hdf_multi1d_cmplx_wrt_F");
  hdf_multi1d_cmplx_wrt_F(pion_PDF_file, lngth, tr_prpgtrs_fft);
  pion_PDF_file.pop();

PASS_ln;
  hdf_multi1d_cmplx_wrt_class cmplx_wrt(pion_PDF_file,lngth);
  cmplx_wrt.wrt("hdf_multi1d_cmplx_wrt_class",tr_prpgtrs_fft);
  cmplx_wrt.wrt("hdf_multi1d_cmplx_wrt_class_1",tr_prpgtrs_fft);

  cmplx_wrt.wrt("hdf_multi1d_cmplx_wrt_class_fft_class",tr_prpgtrs_fft_class);
   
/*
PASS_ln;
  pion_PDF_file.push("hdf_multi1d_cmplx_wrt_class_2");
  cmplx_wrt.wrt("",tr_prpgtrs_fft);
  pion_PDF_file.pop();

PASS_ln;
  pion_PDF_file.push("hdf_multi1d_cmplx_wrt_class_3");
  cmplx_wrt.wrt(tr_prpgtrs_fft);
  pion_PDF_file.pop();
*/

PASS_ln;
  
//  pion_PDF_file.push("hdf_multi1d_cmplx_wrt_F_Class_2");
//  hdf_multi1d_cmplx_wrt_F(pion_PDF_file,lngth,tr_prpgtrs_fft);
//  pion_PDF_file.pop();
PASS_ln;  
  hdf_multi1d_cmplx_wrt_F1(pion_PDF_file,"hdf_multi1d_cmplx_wrt_F1",lngth,tr_prpgtrs_fft);
PASS_ln;
  pion_PDF_file.close();

  cout << "==========================\n"
       << endl;

  cout << "nrow = " << nrow[0] << " " << nrow[1]
       << " " << nrow[2] << " " << nrow[3] << "\n"
       << endl;

  cout << "bndry_cndtn = " << bndry_cndtns[0] << " " << bndry_cndtns[1]
       << " " << bndry_cndtns[2] << " " << bndry_cndtns[3] << "\n"
       << endl;

  cout << pi_mtrx_elmnt << "\n"
       << endl;

  cout << "==========================" << endl;

  END_CODE();
  finalize();
  exit(0);
}
