#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "gauge_link.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl

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

  LatticePropagatorReImPrt quark_prpgtrRe = real(quark_prpgtr);
  LatticePropagatorReImPrt quark_prpgtrIm = imag(quark_prpgtr);

  HDF5Writer quark_prpgtr_1st_file;
  quark_prpgtr_1st_file.open("quark_prpgtr_1st.h5", HDF5Base::ate);
  quark_prpgtr_1st_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
  quark_prpgtr_1st_file.write("Re", quark_prpgtrRe, HDF5Base::trunc);
  quark_prpgtr_1st_file.write("Im", quark_prpgtrIm, HDF5Base::trunc);
  quark_prpgtr_1st_file.close();
//============== 1st propagator sink smearing (not applied here)!!!!================ 
/*  try
  {
    // Get the name of the smearing
    string smearing_name;
    read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Sink/Sink_Smearing/SinkType", smearing_name);

    QDPIO::cout << "sink smearing:" << smearing_name << endl;

    Handle<QuarkSourceSink<LatticePropagator>>
        sinkSmearing(ThePropSinkSmearingFactory::Instance().createObject(smearing_name, xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Sink/Sink_Smearing", u));

    (*sinkSmearing)(quark_prpgtr);
  }
  catch (const string &e)
  {
    QDPIO::cerr << ": Caught Exception : " << e << endl;
    QDP_abort(1);
  }
*/
  //======================Feynman-Hellman Source=========================

  //Src^{H}(Delta)=gaugelink(z,z+\Delta \hat \mu)*quark_prpatr(z+Delta), Needs to shift quark_prpatr(z) to quark_prpatr quark_prpatr(z+Delta)

  //calculate propagator with FH source
  LatticePropagator FH_quark_prpgtr = zero;
  try
  {
    LatticeColorMatrix gaugelink = gauge_link(u, +1, dir, steps);

    LatticePropagator quark_prpgtr_Delta = quark_prpgtr;
    LatticePropagator quark_prpgtr_Delta_tmp = quark_prpgtr;

    if (steps > 0)
    {
      for (int i = 1; i <= steps; i++)
      {
        quark_prpgtr_Delta_tmp = shift(quark_prpgtr_Delta, +1, dir);
        quark_prpgtr_Delta = quark_prpgtr_Delta_tmp;
      }
    }
    else if (steps < 0)
    {
      for (int i = -1; i >= steps; i--)
      {
        quark_prpgtr_Delta_tmp = shift(quark_prpgtr_Delta, -1, dir);
        quark_prpgtr_Delta = quark_prpgtr_Delta_tmp;
      }
    }
    else
    {
      quark_prpgtr_Delta = quark_prpgtr;
    }
    //Gamma(8)=gamma^3 comes from the definition of quasi-PDF

    //  LatticePropagator FH_src = Gamma(8)* gaugelink * quark_prpgtr_Delta; Wrong!, QDP does not support Gamma multiply with ColorMatrix

   LatticePropagator FH_src = gaugelink * (Gamma(8) * quark_prpgtr_Delta);


   LatticePropagatorReImPrt quark_prpgtr_shftd_Re = real(quark_prpgtr_Delta);
   LatticePropagatorReImPrt quark_prpgtr_shftd_Im = imag(quark_prpgtr_Delta);

   HDF5Writer quark_prpgtr_shftd_file;
   quark_prpgtr_shftd_file.open("quark_prpgtr_shftd.h5", HDF5Base::ate);
   quark_prpgtr_shftd_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
   quark_prpgtr_shftd_file.write("Re", quark_prpgtr_shftd_Re, HDF5Base::trunc);
   quark_prpgtr_shftd_file.write("Im", quark_prpgtr_shftd_Im, HDF5Base::trunc);
   quark_prpgtr_shftd_file.close();


  //write out gauge link
    LatticeColorMatrixReImPrt gaugelink_Re = real(gaugelink);
    LatticeColorMatrixReImPrt gaugelink_Im = imag(gaugelink);

    HDF5Writer gaugelink_file;
    gaugelink_file.open("gaugelink.h5", HDF5Base::ate);
    gaugelink_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
    gaugelink_file.write("Re", gaugelink_Re, HDF5Base::trunc);
    gaugelink_file.write("Im", gaugelink_Im, HDF5Base::trunc);
    gaugelink_file.close();

   
 //write out Feynman-Hellman source
    LatticePropagatorReImPrt FH_src_Re = real(FH_src);
    LatticePropagatorReImPrt FH_src_Im = imag(FH_src);

    HDF5Writer FH_src_file;
    FH_src_file.open("FH_src.h5", HDF5Base::ate);
    FH_src_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
    FH_src_file.write("Re", FH_src_Re, HDF5Base::trunc);
    FH_src_file.write("Im", FH_src_Im, HDF5Base::trunc);
    FH_src_file.close();




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

  LatticePropagatorReImPrt FH_quark_prpgtrRe = real(FH_quark_prpgtr);
  LatticePropagatorReImPrt FH_quark_prpgtrIm = imag(FH_quark_prpgtr);

  HDF5Writer FH_quark_prpgtr_file;
  FH_quark_prpgtr_file.open("FH_quark_prpgtr.h5", HDF5Base::ate);
  FH_quark_prpgtr_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
  FH_quark_prpgtr_file.write("Re", FH_quark_prpgtrRe, HDF5Base::trunc);FH_quark_prpgtr_file.write("Im", FH_quark_prpgtrIm, HDF5Base::trunc);
  FH_quark_prpgtr_file.close();

//=========================pion-correlation function========================

  LatticePropagator anti_quark_prpgtr = Gamma(15) * quark_prpgtr * Gamma(15);
  LatticeComplex tr_prpgtrs = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * FH_quark_prpgtr);

  LatticeDouble tr_prpgtrs_Re = real(tr_prpgtrs);
  LatticeDouble tr_prpgtrs_Im = imag(tr_prpgtrs);

  HDF5Writer tr_prpgtrs_file;
  tr_prpgtrs_file.open("tr_prpgtrs.h5", HDF5Base::ate);
  tr_prpgtrs_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
  tr_prpgtrs_file.write("Re", tr_prpgtrs_Re, HDF5Base::trunc);
  tr_prpgtrs_file.write("Im", tr_prpgtrs_Im, HDF5Base::trunc);
  tr_prpgtrs_file.close();


  END_CODE();
  finalize();
  exit(0);
}
