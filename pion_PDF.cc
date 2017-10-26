#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "gauge_link.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl

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

  //read in
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

  XMLReader gauge_file_xml, gauge_xml;
  multi1d<LatticeColorMatrix> u(Nd);

  //HotSt(u);

  gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);
  unitarityCheck(u);

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

  try
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

  LatticeColorMatrix g = one;
  //  rgauge(u,g);

  //calculate propagator with FH source
  LatticePropagator FH_quark_prpgtr = zero;
  try
  {
    LatticeColorMatrix gaugelink = gauge_link(u, +1, dir, steps);

    //quark_prpgtr = g * quark_prpgtr;

    //======================Feynman-Hellman Source=========================

    //Src^{H}(Delta)=gaugelink(z,z+\Delta \hat \mu)*quark_prpatr(z+Delta), Needs to shift quark_prpatr(z) to quark_prpatr quark_prpatr(z+Delta)

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

    //======================Feynman-Hellman Source=========================

    Handle<FermState<T, P, P>> state(S.createState(u));

    LinOpSysSolverCGEnv::registerAll();

    QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

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

  LatticePropagator anti_quark_prpgtr = Gamma(15) * quark_prpgtr * Gamma(15);
  LatticeComplex tr_prpgtrs = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * FH_quark_prpgtr);

  SftMom phases(0, true, Nd - 1);

  multi2d<DComplex> hsum;

  hsum = phases.sft(tr_prpgtrs);

  XMLFileWriter xml_out("pion_correlator");
  XMLArrayWriter momenta(xml_out, phases.numMom());

  push(momenta, "PseudoScalar"); // Array will be called PseudoScalar

  for (int i = 0; i < phases.numMom(); i++)
  {
    push(momenta);
    write(momenta, "mom_index", i);
    write(momenta, "mom", phases.numToMom(i));
    write(momenta, "correlator", hsum[i]);
    pop(momenta);
  }
  pop(momenta);

  xml_out.flush();
  //pop(xml_out);

  END_CODE();
  finalize();
  exit(0);
}
