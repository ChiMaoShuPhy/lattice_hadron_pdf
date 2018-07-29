#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "barspinmat_w.cc"
#include "gauge_link.h"
#include "multi_shift.h"

using namespace std;
using namespace QDP;
using namespace Chroma;
using namespace BaryonSpinMats;

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

  //=============  rgauge  ==============;
//   rgauge(u);

  //======================================================
  //

  //
  //======================================================

  //=================== 1st inversion =====================

  Handle<FermState<T, P, P>> state(S.createState(u));
  LinOpSysSolverCGEnv::registerAll();
  QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

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

  LatticePropagatorReImPrt quark_prpgtrRe = real(quark_prpgtr);
  LatticePropagatorReImPrt quark_prpgtrIm = imag(quark_prpgtr);

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
  LatticeColorMatrix gaugelink = zero;
  LatticePropagator quark_prpgtr_Delta = zero;
  LatticePropagator FH_src = zero;

  try
  {
    gaugelink = gauge_link(u, +1, dir, steps);
    /*
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
*/
    //Gamma(8)=gamma^3 comes from the definition of quasi-PDF

    quark_prpgtr_Delta = multi_shift<LatticePropagator>(quark_prpgtr, steps, dir, nrow[dir], bndry_cndtns[dir]);

    //  LatticePropagator FH_src = Gamma(8)* gaugelink * quark_prpgtr_Delta; Wrong!, QDP does not support Gamma multiply with ColorMatrix

    //========!!!!!!! CAUTION HERE, gamma^z in bi-local operator insertion already multiplied here in the Feynman-Hellman source !!!!!!!==========
    
     FH_src = gaugelink * (Gamma(8) * quark_prpgtr_Delta);

//drop gamma(8) for test: zero dispalcement FH_prpgtr should be identical to quark prpgtr
//    FH_src = gaugelink * ( quark_prpgtr_Delta );

    //======================Feynman-Hellman Source=========================

    Handle<FermState<T, P, P>> state(S.createState(u));

    LinOpSysSolverCGEnv::registerAll();

    QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

    //===============2nd inversion Feynman-Hellman propagator==================

    XMLFileWriter xml_out("fh_xmldat");
    push(xml_out, "FH_proton_crrltr");
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
 
 
  //=========================proton-correlation function========================
  // reference from baryon_3pt.tex in chroma folder
  // parity progector
    SpinMatrix cg5P_pstv = Cg5();
  
  //!!!! IMPORTNT
    SpinMatrix T_unpol = Tunpol();
  //!!!! IMPORTNT


  // sum_x <P(t,y)| u-bar(x)W[x,x+delta_z]u(x+delta_z) |P(0,0)> contraction
  LatticeComplex u_tr_prpgtrs_1 =
      trace(T_unpol * traceColor(FH_quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

  LatticeComplex u_tr_prpgtrs_2 =
      trace(T_unpol * traceColor(FH_quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

  LatticeComplex u_tr_prpgtrs_3 =
      trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * FH_quark_prpgtr))));

  LatticeComplex u_tr_prpgtrs_4 =
      trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * FH_quark_prpgtr)));

  // sum_x <P(t,y)| d-bar(x)W[x,x+delta_z]d(x+delta_z) |P(0,0)> contraction
//!!! REPLACE FH_quark_prpgatr with quark_prpgtr to check whether d quark reduce to 2-pt function

/*  LatticeComplex d_tr_prpgtrs_1 =
      trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

  LatticeComplex d_tr_prpgtrs_2 =
      trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));
*/
  LatticeComplex d_tr_prpgtrs_1 =
      trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(FH_quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

  LatticeComplex d_tr_prpgtrs_2 =
      trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(FH_quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

//=============================================================================  



//proton-prton two point correlation
   LatticeComplex TwoPnt_tr_prpgtrs_1 =
      trace( T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

   LatticeComplex TwoPnt_tr_prpgtrs_2 =
      trace( T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));
  
  LatticeComplex TwoPnt_crrltn_fn = TwoPnt_tr_prpgtrs_1 + TwoPnt_tr_prpgtrs_2;

/*
  SftMom phases(0, true, Nd-1);

  multi2d<DComplex> hsum;
  hsum = phases.sft(TwoPnt_crrltn_fn);


  cout<<"================================="<<endl;
   for(int mom = 0; mom < phases.numMom(); mom++) 
   {
    QDPIO::cout << "mom " << mom << endl;
    for(int t = 0; t < hsum.size1(); t++)
    {
      QDPIO::cout << "t " << t << hsum[mom][t]<<endl;
    }
  }
*/
//take the real, imaginary part

  LatticePropagatorReImPrt quark_prpgtr_Re = real(quark_prpgtr);
  LatticePropagatorReImPrt quark_prpgtr_Im = imag(quark_prpgtr);

  LatticePropagatorReImPrt quark_prpgtr_shftd_Re = real(quark_prpgtr_Delta);
  LatticePropagatorReImPrt quark_prpgtr_shftd_Im = imag(quark_prpgtr_Delta);

  LatticeColorMatrixReImPrt gaugelink_Re = real(gaugelink);
  LatticeColorMatrixReImPrt gaugelink_Im = imag(gaugelink);

  LatticePropagatorReImPrt FH_src_Re = real(FH_src);
  LatticePropagatorReImPrt FH_src_Im = imag(FH_src);

  LatticePropagatorReImPrt FH_quark_prpgtrRe = real(FH_quark_prpgtr);
  LatticePropagatorReImPrt FH_quark_prpgtrIm = imag(FH_quark_prpgtr);

  LatticeDouble u_tr_prpgtrs_Re = real(u_tr_prpgtrs_1 + u_tr_prpgtrs_2 + u_tr_prpgtrs_3 + u_tr_prpgtrs_4);
  LatticeDouble u_tr_prpgtrs_Im = imag(u_tr_prpgtrs_1 + u_tr_prpgtrs_2 + u_tr_prpgtrs_3 + u_tr_prpgtrs_4);

  LatticeDouble d_tr_prpgtrs_Re = real(d_tr_prpgtrs_1 + d_tr_prpgtrs_2);
  LatticeDouble d_tr_prpgtrs_Im = imag(d_tr_prpgtrs_1 + d_tr_prpgtrs_2);

  LatticeDouble TwoPnt_tr_prpgtrs_Re = real(TwoPnt_crrltn_fn);
  LatticeDouble TwoPnt_tr_prpgtrs_Im = imag(TwoPnt_crrltn_fn);
  


  string mtrx_elmnt = "<P(y,t)|q-bar(x).ga_3.W[x,x+" + to_string(steps) + "_stps_in_" + to_string(dir) + "].q(x+" + to_string(steps) + "_stps_in_" + to_string(dir) + ")|P(0,0)>";


  HDF5Writer proton_PDF_file;
  proton_PDF_file.open("proton_PDF.h5", HDF5Base::ate);
  proton_PDF_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.

  proton_PDF_file.push(mtrx_elmnt);
  
  proton_PDF_file.push("quark_prpgtrs");
  proton_PDF_file.write("Re", quark_prpgtr_Re, HDF5Base::trunc);
  proton_PDF_file.write("Im", quark_prpgtr_Im, HDF5Base::trunc);
  proton_PDF_file.pop();

  proton_PDF_file.push("quark_shftd_prpgtrs");
  proton_PDF_file.write("Re", quark_prpgtr_shftd_Re, HDF5Base::trunc);
  proton_PDF_file.write("Im", quark_prpgtr_shftd_Im, HDF5Base::trunc);
  proton_PDF_file.pop();
  

  proton_PDF_file.push("gauge_link");
  proton_PDF_file.write("Re", gaugelink_Re, HDF5Base::trunc);
  proton_PDF_file.write("Im", gaugelink_Im, HDF5Base::trunc);
  proton_PDF_file.pop();

  proton_PDF_file.push("FH_src");
  proton_PDF_file.write("Re", FH_src_Re, HDF5Base::trunc);
  proton_PDF_file.write("Im", FH_src_Im, HDF5Base::trunc);
  proton_PDF_file.pop();
  
  proton_PDF_file.push("FH_quark_prpgtr");
    proton_PDF_file.write("Re", FH_quark_prpgtrRe, HDF5Base::trunc);
    proton_PDF_file.write("Im", FH_quark_prpgtrIm, HDF5Base::trunc);
  proton_PDF_file.pop();


  proton_PDF_file.push("u_tr_prpgtrs");
    proton_PDF_file.write("Re", u_tr_prpgtrs_Re, HDF5Base::trunc);
    proton_PDF_file.write("Im", u_tr_prpgtrs_Im, HDF5Base::trunc);
  proton_PDF_file.pop();

  proton_PDF_file.push("d_tr_prpgtrs");
    proton_PDF_file.write("Re", d_tr_prpgtrs_Re, HDF5Base::trunc);
    proton_PDF_file.write("Im", d_tr_prpgtrs_Im, HDF5Base::trunc);
  proton_PDF_file.pop();

  proton_PDF_file.push("TwoPnt_tr_prpgtrs");
    proton_PDF_file.write("Re", TwoPnt_tr_prpgtrs_Re, HDF5Base::trunc);
    proton_PDF_file.write("Im", TwoPnt_tr_prpgtrs_Im, HDF5Base::trunc);
  proton_PDF_file.pop();

  proton_PDF_file.pop();
  
  proton_PDF_file.close();
  
  cout<<"==========================\n"<<endl;

  cout<<"nrow = "<<nrow[0]<<" "<<nrow[1]
      <<" "<<nrow[2]<<" "<<nrow[3]<<"\n"<<endl;

  cout<<"bndry_cndtn = "<<bndry_cndtns[0]<<" "<<bndry_cndtns[1]
      <<" "<<bndry_cndtns[2]<<" "<<bndry_cndtns[3]<<"\n"<<endl;

  cout<< mtrx_elmnt <<"\n"<<endl;

  cout<<"=========================="<<endl;

  END_CODE();
  finalize();
  exit(0);
}
