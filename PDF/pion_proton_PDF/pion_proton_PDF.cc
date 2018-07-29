#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "barspinmat_w.cc"
#include "gauge_link.h"
#include "multi_shift.h"
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
    pion_file.open("pion_PDF.h5", HDF5Base::ate);
    pion_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.

    HDF5Writer proton_file;
    proton_file.open("proton_PDF.h5", HDF5Base::ate);
    proton_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.
//--1push
    pion_file.push(cfg_file_name);
    proton_file.push(cfg_file_name);

  for (int src_iter = 0; src_iter < sobol_src.size(); src_iter++)
  {
    for (int i = 0; i < Nd; i++)
    {
      t_src[i] = sobol_src[src_iter].int_arry[i];
    }
//============= 1st inversion ===============

  Handle<FermState<T, P, P>> state(S.createState(u));
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


    //================Feynman-Hellman Source============

    //Src^{H}(Delta)=gaugelink(z,z+\Delta \hat \mu)*quark_prpatr(z+Delta), Needs to shift quark_prpatr(z) to quark_prpatr quark_prpatr(z+Delta)

    //calculate propagator with FH source
    LatticePropagator FH_quark_prpgtr = zero;
    LatticeColorMatrix gaugelink = zero;
    LatticePropagator quark_prpgtr_Delta = zero;
    LatticePropagator FH_src = zero;
//Sobol sequence as source position
    string src_pstn = to_string(sobol_src[src_iter].int_arry[0]) + " " + to_string(sobol_src[src_iter].int_arry[1]) + " " + to_string(sobol_src[src_iter].int_arry[2]) + " " + to_string(sobol_src[src_iter].int_arry[3]) + " ";
//--2push    
    pion_file.push(src_pstn);
    proton_file.push(src_pstn);
//   hdf_multi1d_cmplx_wrt pion_wrt(pion_file,t_lngth);
//   hdf_multi1d_cmplx_wrt //proton_wrt(//proton_file,t_lngth);

    multi1d<int> hdrn_mmntm(Nd);
    for (int dir_mu = 0; dir_mu < Nd; dir_mu++)
    {
      hdrn_mmntm[dir_mu] = 0;
    }

    //loop over all sptial direction
    //mu <P_mu|psi(r+x_mu)*W(r+x_mu,r)*gamma_mu*psi(r)|P_mu>
//    for (int dir_mu = 0; dir_mu < 1; dir_mu++)  
    for (int dir_mu = 0; dir_mu < Nd-1; dir_mu++)
    {
      if (dir_mu != j_decay)
      {
        //loop dispalcement step in PDF bilocal operator
//        for (int steps = 0; steps < 2; steps++)
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
          if (steps == 0)//only need to calculate two-point function once!!!!
          {
            //pion-pion two point correlation
            LatticeComplex pion_TwoPnt_tr_prpgtrs = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * quark_prpgtr);

            //proton-prton two point correlation
            LatticeComplex proton_TwoPnt_tr_prpgtrs_1 =
                trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

            LatticeComplex proton_TwoPnt_tr_prpgtrs_2 =
                trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

            FFT pion_TwoPnt_tr_prpgtrs_fft(pion_TwoPnt_tr_prpgtrs, j_decay, t_lngth);
            FFT proton_TwoPnt_tr_prpgtrs_1_fft(proton_TwoPnt_tr_prpgtrs_1, j_decay, t_lngth);
            FFT proton_TwoPnt_tr_prpgtrs_2_fft(proton_TwoPnt_tr_prpgtrs_2, j_decay, t_lngth);

            multi1d<DComplex> pion_TwoPnt_tr_prpgtrs_fftd;
            multi1d<DComplex> proton_TwoPnt_tr_prpgtrs_1_fftd;
            multi1d<DComplex> proton_TwoPnt_tr_prpgtrs_2_fftd;           

            multi1d<int> hdrn_mmntm_zero(4);
            hdrn_mmntm_zero[0] = 0; hdrn_mmntm_zero[1] = 0; 
            hdrn_mmntm_zero[2] = 0; hdrn_mmntm_zero[3] = 0; 

            pion_TwoPnt_tr_prpgtrs_fftd = pion_TwoPnt_tr_prpgtrs_fft.prjctn(hdrn_mmntm_zero);

            proton_TwoPnt_tr_prpgtrs_1_fftd = proton_TwoPnt_tr_prpgtrs_1_fft.prjctn(hdrn_mmntm_zero);
            proton_TwoPnt_tr_prpgtrs_2_fftd = proton_TwoPnt_tr_prpgtrs_2_fft.prjctn(hdrn_mmntm_zero);
    
            hdf_multi1d_cmplx_wrt_F(pion_file,"TwoPnt_tr_prpgtrs_fft",t_lngth,pion_TwoPnt_tr_prpgtrs_fftd);
             
//--3push  
            proton_file.push("TwoPnt_tr_prpgtrs");      
            
            hdf_multi1d_cmplx_wrt_F(proton_file,"TwoPnt_tr_prpgtrs_1_fft",t_lngth,proton_TwoPnt_tr_prpgtrs_1_fftd);
            
            hdf_multi1d_cmplx_wrt_F(proton_file,"TwoPnt_tr_prpgtrs_2_fft",t_lngth,proton_TwoPnt_tr_prpgtrs_2_fftd);
//-3pop               
            proton_file.pop();
          }


          // sum_x <P(t,y)| u-bar(x)W[x,x+delta_z]u(x+delta_z) |P(0,0)> contraction
          // pion trced propagator
          LatticeComplex pion_tr_prpgtrs = trace(Gamma(15) * adj(anti_quark_prpgtr) * Gamma(15) * FH_quark_prpgtr);


          // proton trced propagator
          LatticeComplex proton_u_tr_prpgtrs_1 =
              trace(T_unpol * traceColor(FH_quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

          LatticeComplex proton_u_tr_prpgtrs_2 =
              trace(T_unpol * traceColor(FH_quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

          LatticeComplex proton_u_tr_prpgtrs_3 =
              trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * FH_quark_prpgtr))));

          LatticeComplex proton_u_tr_prpgtrs_4 =
              trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(quark_prpgtr * cg5P_pstv, cg5P_pstv * FH_quark_prpgtr)));

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
//Fourier Transform, overall factor 1/sqrt(L_1*L_2*L3) contained in Fourier              multi1d<DComplex> pion_tr_prpgtrs_fft_all(t_lngth);
          FFT pion_tr_prpgtrs_fft(pion_tr_prpgtrs, j_decay, t_lngth);

          FFT proton_u_tr_prpgtrs_1_fft(proton_u_tr_prpgtrs_1, j_decay, t_lngth);
          FFT proton_u_tr_prpgtrs_2_fft(proton_u_tr_prpgtrs_2, j_decay, t_lngth);
          FFT proton_u_tr_prpgtrs_3_fft(proton_u_tr_prpgtrs_3, j_decay, t_lngth);
          FFT proton_u_tr_prpgtrs_4_fft(proton_u_tr_prpgtrs_4, j_decay, t_lngth);

          FFT proton_d_tr_prpgtrs_1_fft(proton_d_tr_prpgtrs_1, j_decay, t_lngth);
          FFT proton_d_tr_prpgtrs_2_fft(proton_d_tr_prpgtrs_2, j_decay, t_lngth);
         
          multi1d<DComplex> pion_tr_prpgtrs_fftd(t_lngth);

          multi1d<DComplex> proton_u_tr_prpgtrs_1_fftd(t_lngth);
          multi1d<DComplex> proton_u_tr_prpgtrs_2_fftd(t_lngth);
          multi1d<DComplex> proton_u_tr_prpgtrs_3_fftd(t_lngth);
          multi1d<DComplex> proton_u_tr_prpgtrs_4_fftd(t_lngth);

          multi1d<DComplex> proton_d_tr_prpgtrs_1_fftd(t_lngth);
          multi1d<DComplex> proton_d_tr_prpgtrs_2_fftd(t_lngth);            

          string mtrx_elmnt = "<P(y,t)|qbr(x).ga_" + to_string(dir_mu + 1) + ".W" + ".q(x+" + to_string(steps) + "_stps_in_" + to_string(dir_mu) + ")|P(0,0)>";
//--4push
          pion_file.push(mtrx_elmnt);
          proton_file.push(mtrx_elmnt);

//           for (int mom = 1; mom < 2; mom++)
          for (int mom = 0; mom < nrow[dir_mu]; mom++)
          {
            hdrn_mmntm[0] = 0; hdrn_mmntm[1] = 0; hdrn_mmntm[2] = 0; hdrn_mmntm[3] = 0; 

            hdrn_mmntm[dir_mu] = mom;

            pion_tr_prpgtrs_fftd = pion_tr_prpgtrs_fft.prjctn(hdrn_mmntm);
            
            proton_u_tr_prpgtrs_1_fftd = proton_u_tr_prpgtrs_1_fft.prjctn(hdrn_mmntm);
            proton_u_tr_prpgtrs_2_fftd = proton_u_tr_prpgtrs_2_fft.prjctn(hdrn_mmntm);
            proton_u_tr_prpgtrs_3_fftd = proton_u_tr_prpgtrs_3_fft.prjctn(hdrn_mmntm);
            proton_u_tr_prpgtrs_4_fftd = proton_u_tr_prpgtrs_4_fft.prjctn(hdrn_mmntm);

            proton_d_tr_prpgtrs_1_fftd = proton_d_tr_prpgtrs_1_fft.prjctn(hdrn_mmntm);
            proton_d_tr_prpgtrs_2_fftd = proton_d_tr_prpgtrs_2_fft.prjctn(hdrn_mmntm);

            string mom_string = "{" + to_string(hdrn_mmntm[0]) + "," + to_string(hdrn_mmntm[1]) + "," + to_string(hdrn_mmntm[2]) + "," + to_string(hdrn_mmntm[3]) + "}";

 //--5push           
            pion_file.push(mom_string);            
            proton_file.push(mom_string);
           
             hdf_multi1d_cmplx_wrt_F(pion_file,"tr_prpgtrs_fft",t_lngth,pion_tr_prpgtrs_fftd);
 
//           pion_wrt.wrt("tr_prpgtrs_fft",pion_tr_prpgtrs_fftd);

 //--6push
            proton_file.push("u_tr_prpgtrs_fft");
            hdf_multi1d_cmplx_wrt_F(proton_file,"1st",t_lngth,proton_u_tr_prpgtrs_1_fftd);
            hdf_multi1d_cmplx_wrt_F(proton_file,"2nd",t_lngth,proton_u_tr_prpgtrs_2_fftd);
            hdf_multi1d_cmplx_wrt_F(proton_file,"3rd",t_lngth,proton_u_tr_prpgtrs_3_fftd);
            hdf_multi1d_cmplx_wrt_F(proton_file,"4th",t_lngth,proton_u_tr_prpgtrs_4_fftd);
 //--6pop            
            proton_file.pop();
 //--7push
            proton_file.push("d_tr_prpgtrs_fft");
            hdf_multi1d_cmplx_wrt_F(proton_file,"1st",t_lngth,proton_d_tr_prpgtrs_1_fftd);
            hdf_multi1d_cmplx_wrt_F(proton_file,"2nd",t_lngth,proton_d_tr_prpgtrs_2_fftd);
 //--7pop            
            proton_file.pop();

 //--5pop            
            pion_file.pop();//pop mom_string
            
            proton_file.pop();//pop mom_string
          }
 //--4pop        
          pion_file.pop();// pop mtrx_elmnt
          proton_file.pop();// pop mtrx_elmnt

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
    }
//--2pop
    pion_file.pop();   //pop source position
    proton_file.pop(); //pop source position
  }
//--1pop
  pion_file.pop();   //pop cfg file name
  proton_file.pop(); //pop cfg file name

  pion_file.close();
  proton_file.close();
  END_CODE();
  finalize();
  exit(0);
}