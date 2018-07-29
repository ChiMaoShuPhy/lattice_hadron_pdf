#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "barspinmat_w.cc"
#include "gauge_link.h"
#include "multi_shift.h"

#include "timeslice_FT.h"

#include "sobol.h"
#include "sobol.cpp"
#include <vector>

//with sobol sequence as source position

using namespace std;
using namespace QDP;
using namespace Chroma;
using namespace BaryonSpinMats;

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl

typedef OLattice<PSpinMatrix<PColorMatrix<RScalar<double>, Nc>, Ns>> LatticePropagatorReImPrt;

typedef OLattice<PScalar<PColorMatrix<RScalar<double>, Nc>>> LatticeColorMatrixReImPrt;

struct int4
{
  int int_arry[4];
};

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

void sobol_src_gnrtr(int srnd_seed, int lngth, multi1d<int> nrow,
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

  for (unsigned long long i = 0; i < lngth + seed; ++i)
  {
    // Print a few dimensions of each point.

    for (unsigned d = 0; d < 4; ++d)
    {
      const int s = (int)(nrow[d] * sobol::sample(i, d));
      if (skip < seed)
      {
        ;
      }
      else
      {
        src_tmp.int_arry[d] = s;
      }
    }
    if (skip < seed)
    {
      ;
    }
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
  // multi1d<int> t_srce;

  //  int dir;
  //  int steps;

  Real wvf_param;
  int wvfIntPar;

  //================= set parameters ===============

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
    // read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Source_Smearing/t_srce", t_srce);
    //    read(xml_in, "/chroma/Param/InlineMeasurements/elem/displacement/dir", dir);
    //  read(xml_in, "/chroma/Param/InlineMeasurements/elem/displacement/steps", steps);
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

  //============================
  //

  //
  //============================

  //============= 1st inversion ===============

  Handle<FermState<T, P, P>> state(S.createState(u));
  LinOpSysSolverCGEnv::registerAll();
  QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

  multi1d<int> t_src(Nd);

  int probe = 0;

  // generating sobol sequence as source position//

  vector<int4> sobol_src;
  sobol_src_gnrtr(Sobol_src_seed, Sobol_src_Lngth, nrow, sobol_src);

  cout << "========================" << endl;
  for (int iter = 0; iter < sobol_src.size(); iter++)
  {
    cout << "!!!" << iter << "!!!" << endl;
    cout << sobol_src[iter].int_arry[0] << " " << sobol_src[iter].int_arry[1] << " "
         << sobol_src[iter].int_arry[2] << " " << sobol_src[iter].int_arry[3] << endl;
  }

  for (int src_iter = 0; src_iter < sobol_src.size(); src_iter++)
  {
    for (int i = 0; i < Nd; i++)
    {
      t_src[i] = sobol_src[src_iter].int_arry[i];
    }

    //    cout<<"*******"<<t_src[0]<<t_src[1]<<t_src[2]<<t_src[3]<<endl;

    quark_src = getPointSource(t_src, j_decay);

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

    //================Feynman-Hellman Source============

    //Src^{H}(Delta)=gaugelink(z,z+\Delta \hat \mu)*quark_prpatr(z+Delta), Needs to shift quark_prpatr(z) to quark_prpatr quark_prpatr(z+Delta)

    //calculate propagator with FH source
    LatticePropagator FH_quark_prpgtr = zero;
    LatticeColorMatrix gaugelink = zero;
    LatticePropagator quark_prpgtr_Delta = zero;
    LatticePropagator FH_src = zero;

    HDF5Writer pion_PDF_file;
    pion_PDF_file.open("pion_PDF.h5", HDF5Base::ate);
    pion_PDF_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.

    HDF5Writer proton_PDF_file;
    proton_PDF_file.open("proton_PDF.h5", HDF5Base::ate);
    proton_PDF_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.

    pion_PDF_file.push(cfg_file_name);
    proton_PDF_file.push(cfg_file_name);

    string src_pstn = to_string(sobol_src[src_iter].int_arry[0]) + " " + to_string(sobol_src[src_iter].int_arry[1]) + " " + to_string(sobol_src[src_iter].int_arry[2]) + " " + to_string(sobol_src[src_iter].int_arry[3]) + " ";
    pion_PDF_file.push(src_pstn);
    proton_PDF_file.push(src_pstn);
          
    multi1d<int> hdrn_mmntm(Nd);
    for (int dir_mu = 0; dir_mu < Nd; dir_mu++)
    {
      hdrn_mmntm[dir_mu] = 0;
    }
          
    // gamma matrices

    //loop over all sptial direction
    //mu <P_mu|psi(r+x_mu)*W(r+x_mu,r)*gamma_mu*psi(r)|P_mu>
    for (int dir_mu = 0; dir_mu < Nd; dir_mu++)
    {
      if (dir_mu != j_decay)
      {
        //loop dispalcement step in PDF bilocal operator
        for (int steps = 0; steps < nrow[dir_mu]; steps++)
        {
          try
          {
            gaugelink = gauge_link(u, +1, dir_mu, steps);

            quark_prpgtr_Delta = multi_shift<LatticePropagator>(quark_prpgtr, steps, dir_mu, nrow[dir_mu], bndry_cndtns[dir_mu]);

            //  LatticePropagator FH_src = Gamma(mu)* gaugelink * quark_prpgtr_Delta; Wrong!, QDP does not support Gamma multiply with ColorMatrix
            // mu MUST MATCH the momentum's direction!!!!
            //eg. gamma_3 ~ P = (0,0,P3,P4), gamma_1 ~ P = (P1,0,0,P4), to increase the statistic.

            //========!!!!!!! CAUTION HERE, gamma_mu in bi-local operator insertion already multiplied here in the Feynman-Hellman source !!!!!!!==========

            FH_src = gaugelink * (Gamma(pow(2, dir_mu + 1)) * quark_prpgtr_Delta);

            //============Feynman-Hellman Source============

            Handle<FermState<T, P, P>> state(S.createState(u));

            LinOpSysSolverCGEnv::registerAll();

            QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

            //===============2nd inversion Feynman-Hellman propagator============

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
          
          swatch.stop();
          //============pion, proton-correlation function==================
          // reference from baryon_3pt.tex in chroma folder
          // parity progector
          SpinMatrix cg5P_pstv = Cg5();

          //!!!! IMPORTNT
          SpinMatrix T_unpol = Tunpol();
          //!!!! IMPORTNT

          LatticePropagator anti_quark_prpgtr = Gamma(15) * quark_prpgtr * Gamma(15);
          
          //loop over dispacement in PDF-bilocal operator
          if (steps == 0)
          //only need to calculate two-point function once!!!!
          {
            //pion-pion two point correlation
            LatticeComplex pion_TwoPnt_tr_prpgtrs = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * quark_prpgtr);

            LatticeDouble pion_TwoPnt_tr_prpgtrs_Re = real(pion_TwoPnt_tr_prpgtrs);
            LatticeDouble pion_TwoPnt_tr_prpgtrs_Im = imag(pion_TwoPnt_tr_prpgtrs);

            //proton-prton two point correlation
            LatticeComplex proton_TwoPnt_tr_prpgtrs_1 =
                trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

            LatticeComplex proton_TwoPnt_tr_prpgtrs_2 =
                trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

            LatticeComplex proton_TwoPnt_crrltn_fn = proton_TwoPnt_tr_prpgtrs_1 + proton_TwoPnt_tr_prpgtrs_2;

            LatticeDouble proton_TwoPnt_tr_prpgtrs_Re = real(proton_TwoPnt_crrltn_fn);
            LatticeDouble proton_TwoPnt_tr_prpgtrs_Im = imag(proton_TwoPnt_crrltn_fn);

            pion_PDF_file.push("TwoPnt_tr_prpgtrs");
            pion_PDF_file.write("Re", pion_TwoPnt_tr_prpgtrs_Re, HDF5Base::trunc);
            pion_PDF_file.write("Im", pion_TwoPnt_tr_prpgtrs_Im, HDF5Base::trunc);
            pion_PDF_file.pop();

            proton_PDF_file.push("TwoPnt_tr_prpgtrs");
            proton_PDF_file.write("Re", proton_TwoPnt_tr_prpgtrs_Re, HDF5Base::trunc);
            proton_PDF_file.write("Im", proton_TwoPnt_tr_prpgtrs_Im, HDF5Base::trunc);
            proton_PDF_file.pop();
          }
          else  ;
                    
          // sum_x <P(t,y)| u-bar(x)W[x,x+delta_z]u(x+delta_z) |P(0,0)> contraction
          // pion trced propagator
          LatticeComplex pion_tr_prpgtrs = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * FH_quark_prpgtr);

          LatticeDouble pion_tr_prpgtrs_Re = real(pion_tr_prpgtrs);
          LatticeDouble pion_tr_prpgtrs_Im = imag(pion_tr_prpgtrs);

          // proton trced propagator
          LatticeComplex proton_u_tr_prpgtrs_1 =
              trace(T_unpol * traceColor(FH_quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

          LatticeComplex proton_u_tr_prpgtrs_2 =
              trace(T_unpol * traceColor(FH_quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

          LatticeComplex proton_u_tr_prpgtrs_3 =
              trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * FH_quark_prpgtr))));

          LatticeComplex proton_u_tr_prpgtrs_4 =
              trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * FH_quark_prpgtr)));
          
          // sum_x <P(t,y)| d-bar(x)W[x,x+delta_z]d(x+delta_z) |P(0,0)> contraction
          //!!! REPLACE FH_quark_prpgatr with quark_prpgtr to check whether d quark reduce to 2-pt function

          /*  LatticeComplex d_tr_prpgtrs_1 =
          trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

          LatticeComplex d_tr_prpgtrs_2 =
           trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr *   cg5P_pstv, cg5P_pstv * quark_prpgtr)));
          */
          LatticeComplex proton_d_tr_prpgtrs_1 =
              trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(FH_quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

          LatticeComplex proton_d_tr_prpgtrs_2 =
              trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(FH_quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

          QDPIO::cout << "\n"
                      << "==============="
                      << "\n"
                      << "Cfg file: " << cfg_file_name << "\n"
                      << "Source position: " << t_src[0] << ", " << t_src[1] << ", " << t_src[2] << ", " << t_src[3] << "\n"
                      << "Correlation function calculation done, time= " << swatch.getTimeInSeconds() << " secs"
                      << "\n"
                      << "===============" << endl;

          LatticeDouble proton_u_tr_prpgtrs_Re = real(proton_u_tr_prpgtrs_1 + proton_u_tr_prpgtrs_2 + proton_u_tr_prpgtrs_3 + proton_u_tr_prpgtrs_4);
          LatticeDouble proton_u_tr_prpgtrs_Im = imag(proton_u_tr_prpgtrs_1 + proton_u_tr_prpgtrs_2 + proton_u_tr_prpgtrs_3 + proton_u_tr_prpgtrs_4);

          LatticeDouble proton_d_tr_prpgtrs_Re = real(proton_d_tr_prpgtrs_1 + proton_d_tr_prpgtrs_2);
          LatticeDouble proton_d_tr_prpgtrs_Im = imag(proton_d_tr_prpgtrs_1 + proton_d_tr_prpgtrs_2);

          string mtrx_elmnt = "<P(y,t)|q-bar(x).ga_" + to_string(dir_mu + 1) + ".W[x,x+" + to_string(steps) + "_stps_in_" + to_string(dir_mu) + "].q(x+" + to_string(steps) + "_stps_in_" + to_string(dir_mu) + ")|P(0,0)>";

          int lngth = nrow[j_decay];

          multi1d<DComplex> pion_tr_prpgtrs_FT_Sum(lngth);
          multi1d<DComplex> proton_u_tr_prpgtrs_FT_Sum(lngth);
          multi1d<DComplex> proton_d_tr_prpgtrs_FT_Sum(lngth);

          multi1d<Real> pion_tr_prpgtrs_FT_Sum_Re(lngth);
          multi1d<Real> pion_tr_prpgtrs_FT_Sum_Im(lngth);
          multi1d<Real> proton_u_tr_prpgtrs_FT_Sum_Re(lngth);
          multi1d<Real> proton_u_tr_prpgtrs_FT_Sum_Im(lngth);
          multi1d<Real> proton_d_tr_prpgtrs_FT_Sum_Re(lngth);
          multi1d<Real> proton_d_tr_prpgtrs_FT_Sum_Im(lngth);

          pion_PDF_file.push(mtrx_elmnt);
          proton_PDF_file.push(mtrx_elmnt);


          for (int mom = 0; mom < nrow[dir_mu]; mom++)
          {
            for (int mu = 0; mu < Nd; mu++)
            {
              hdrn_mmntm[mu] = 0;
            }

            hdrn_mmntm[dir_mu] = mom;

            pion_tr_prpgtrs_FT_Sum =
                timesilce_FT<multi1d<DComplex>, LatticeComplex>(pion_tr_prpgtrs, nrow, hdrn_mmntm, j_decay, Nd);

            proton_u_tr_prpgtrs_FT_Sum =
                timesilce_FT<multi1d<DComplex>, LatticeComplex>(proton_u_tr_prpgtrs_1 + proton_u_tr_prpgtrs_2 + proton_u_tr_prpgtrs_3 + proton_u_tr_prpgtrs_4, nrow, hdrn_mmntm, j_decay, Nd);

            proton_d_tr_prpgtrs_FT_Sum =
                timesilce_FT<multi1d<DComplex>, LatticeComplex>(proton_d_tr_prpgtrs_1 + proton_d_tr_prpgtrs_2,nrow, hdrn_mmntm, j_decay, Nd);
          
            for (int i = 0; i < lngth; i++)
            {

              pion_tr_prpgtrs_FT_Sum_Re[i] = real(pion_tr_prpgtrs_FT_Sum[i]);
              pion_tr_prpgtrs_FT_Sum_Im[i] = imag(pion_tr_prpgtrs_FT_Sum[i]);

              proton_u_tr_prpgtrs_FT_Sum_Re[i] = real(proton_u_tr_prpgtrs_FT_Sum[i]);
              proton_u_tr_prpgtrs_FT_Sum_Im[i] = imag(proton_u_tr_prpgtrs_FT_Sum[i]);

              proton_d_tr_prpgtrs_FT_Sum_Re[i] = real(proton_d_tr_prpgtrs_FT_Sum[i]);
              proton_d_tr_prpgtrs_FT_Sum_Im[i] = imag(proton_d_tr_prpgtrs_FT_Sum[i]);
            }
                       
            string mom_string = "{" + to_string(hdrn_mmntm[0]) + "," + to_string(hdrn_mmntm[1]) + "," + to_string(hdrn_mmntm[2]) + "," + to_string(hdrn_mmntm[3]) + "}";
 
            pion_PDF_file.push(mom_string);
            pion_PDF_file.push("tr_prpgtrs_FT_Sum");
            pion_PDF_file.write("Re", pion_tr_prpgtrs_FT_Sum_Re, HDF5Base::trunc);
            pion_PDF_file.write("Im", pion_tr_prpgtrs_FT_Sum_Im, HDF5Base::trunc);
            pion_PDF_file.pop();
            pion_PDF_file.pop();

            proton_PDF_file.push(mom_string);
            proton_PDF_file.push("u_tr_prpgtrs_FT_Sum");
            proton_PDF_file.write("Re", proton_u_tr_prpgtrs_FT_Sum_Re, HDF5Base::trunc);
            proton_PDF_file.write("Im", proton_u_tr_prpgtrs_FT_Sum_Im, HDF5Base::trunc);
            proton_PDF_file.pop();
            proton_PDF_file.push("d_tr_prpgtrs_FT_Sum");
            proton_PDF_file.write("Re", proton_d_tr_prpgtrs_FT_Sum_Re, HDF5Base::trunc);
            proton_PDF_file.write("Im", proton_d_tr_prpgtrs_FT_Sum_Im, HDF5Base::trunc);
            proton_PDF_file.pop();
            proton_PDF_file.pop();
          }
                    
          //========= Write pion traced propagators =========

          pion_PDF_file.push("tr_prpgtrs");
          pion_PDF_file.write("Re", pion_tr_prpgtrs_Re, HDF5Base::trunc);
          pion_PDF_file.write("Im", pion_tr_prpgtrs_Im, HDF5Base::trunc);
          pion_PDF_file.pop();

          pion_PDF_file.pop();
                  
          //========= Write proton traced propagators =========

          proton_PDF_file.push("u_tr_prpgtrs");
          proton_PDF_file.write("Re", proton_u_tr_prpgtrs_Re, HDF5Base::trunc);
          proton_PDF_file.write("Im", proton_u_tr_prpgtrs_Im, HDF5Base::trunc);
          proton_PDF_file.pop();
          proton_PDF_file.push("d_tr_prpgtrs");
          proton_PDF_file.write("Re", proton_d_tr_prpgtrs_Re, HDF5Base::trunc);
          proton_PDF_file.write("Im", proton_d_tr_prpgtrs_Im, HDF5Base::trunc);
          proton_PDF_file.pop();

          proton_PDF_file.pop();

                    

          cout << "===============\n"

               << "nrow = " << nrow[0] << " " << nrow[1]
               << " " << nrow[2] << " " << nrow[3] << "\n"

               << "bndry_cndtn = " << bndry_cndtns[0] << " " << bndry_cndtns[1]
               << " " << bndry_cndtns[2] << " " << bndry_cndtns[3] << "\n"

               << "Cfg file: " << cfg_file_name << "\n"
               << "Source position: " << t_src[0] << ", " << t_src[1] << ", " << t_src[2] << ", " << t_src[3] << "\n"

               << mtrx_elmnt << "\n"

               << "===============" << endl;
        }
      }
      else
        ;
    }
    pion_PDF_file.pop();   //pop source position
    proton_PDF_file.pop(); //pop source position

    pion_PDF_file.pop();   //pop cfg file name
    proton_PDF_file.pop(); //pop cfg file name

    pion_PDF_file.close();
    proton_PDF_file.close();
  }

  END_CODE();
  finalize();
  exit(0);
}