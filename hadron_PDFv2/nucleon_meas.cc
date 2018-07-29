#include "nucleon_meas.h"
#include "chroma.h"
#include "spin_matrix.h"
#include "barspinmat_w.cc"
#include <iostream>
#include <string>

#include "gauge_link.h"
#include "multi_shift.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl

using namespace QDP;
using namespace std;
using namespace BaryonSpinMats;

typedef OLattice<PSpinMatrix<PColorMatrix<RScalar<double>, Nc>, Ns>> LatticePropagatorReImPrt;

namespace Chroma 
{ 
  // Fill in the blanks
  namespace NucleonEnv 
  { 

    // Function to register with a factory
    // This is boilerplate stuff
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      // Create a Params from the XML
      // Use it to create a Nucleon
      return new Nucleon( Params(xml_in, path) );
    }

    // The name of my measurement for the XML file
    // Change this for each measurement
    const std::string name = "NUCLEON_SYSTEM";

    // Register the measurement in the measurement factory
    namespace { 
      bool registered = false;
    }

    bool registerAll()
    {
      bool success = true;
      if (! registered) { 
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }

    // Param stuff
    Params::Params() { nucleon_system="";frequency = 0; gauge_id=""; smeared_prop_id=""; xml_file=""; filedir=""; consistencyTests = false;}
    
    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
	{
	  XMLReader paramtop(xml_in, path);

	  if (paramtop.count("Frequency") == 1)
	    read(paramtop, "Frequency", frequency);
	  else
	    frequency = 1;
      
	  /* Name of system to do contractions */
	  read(paramtop, "Name", nucleon_system);

	  // Read in the file directory
	  read(paramtop,"Param/file_dir",filedir);

	  // Read in the output propagator/source configuration info
	  read(paramtop, "NamedObject/smeared_prop_id", smeared_prop_id);

	  // Get either the gauge_id in the XML or the ID of the default
	  // field of no explicit ID exists.
	  read(paramtop, "NamedObject/gauge_id", gauge_id);

	  if( paramtop.count("xml_file") != 0 ) { 
	    read(paramtop, "xml_file", xml_file);
	  }
	  else { 
	    xml_file == "";
	  }
	
	}
      catch(const std::string& e) 
	{
	  QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	  QDP_abort(1);
	}

      try {
	XMLReader paramtop(xml_in, path);
      
	// Read in the file directory
	read(paramtop,"Param/ConsistencyTests",consistencyTests);
	if(consistencyTests){
	  QDPIO::cout << "Performing consistency tests" << endl;
	} else {
	  QDPIO::cout << "Not performing consistency tests" << endl;
	};
      }catch(const std::string& e) {
	QDPIO::cout << "No consistency-tests option stated in " << e << " : will NOT perform any consistency tests" << endl;
      }

      // Read your own parameters  from Param/ here
    }


    void
    Params::write(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Write out the file directory
      QDP::write(xml_out,"file_dir",filedir);

      // Write out the output propagator/source configuration info
      //      QDP::write(xml_out, "smeared_prop_id", smeared_prop_id);

      QDP::write(xml_out, "gauge_id", gauge_id);

      if( xml_file != "" ){ 
	QDP::write(xml_out, "xml_file", xml_file);
      }


      pop(xml_out);
    }
  }

  void 
  Nucleon::operator()(long unsigned int update_no,
			  XMLWriter& xml_out) 
  {

    // This bit merely supports providing an external xml_file 
    // for this measurement
    if ( params.xml_file == "" ) { 
      
      func( update_no, xml_out );
    }
    else { 

      // Hash an XML file name from the user specified string
      // and the update number
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      // IN the global output, make a note that the output went
      // to this separate XML file
      push(xml_out, "my_meas");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
  }

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

  void 
  Nucleon::func(unsigned long int update_no, XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch measure_time;
    measure_time.reset();
    measure_time.start();

    PASS(195);

    // Test that the gauge configuration and the propagator we need
    // exist in the map.
    XMLBufferWriter gauge_xml_BW;
    XMLBufferWriter smeared_prop_xml;
    try
      {
	// Try and get at the gauge field if it doesn't exist 
	// an exception will be thrown.
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);
	TheNamedObjMap::Instance().get(params.gauge_id).getRecordXML(gauge_xml_BW);

	// Do the same with the propagator 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.smeared_prop_id);

	// Get its record XML
	TheNamedObjMap::Instance().get(params.smeared_prop_id).getRecordXML(smeared_prop_xml);

      }
    catch( std::bad_cast ) 
      {

	// If the object exists and is of the wrong type we will 
	// land in this catch.
	QDPIO::cerr << params.nucleon_system << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
    catch (const string& e) 
      {
	// If the object is not in the map we will land in this 
	// catch
	QDPIO::cerr << params.nucleon_system << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }

    // If we got here, then the gauge field and prop are both in 
    // the map. Their XML will have been captured.
    // Let us bind the references to a local name so 
    // we don't have to type the long lookup string every time.
    //
    // Note const here means we can't change the field
    //    const multi1d<LatticeColorMatrix>& u = 
    //      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);

    PASS(241);
    
    const LatticePropagator& smeared_prop = 
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.smeared_prop_id);

    XMLReader input_smeared_prop_xml;
    TheNamedObjMap::Instance().get(params.smeared_prop_id).getRecordXML(input_smeared_prop_xml);

    multi1d<int> srce_pt;
    read(input_smeared_prop_xml,"/SinkSmear/PropSource/Source/t_srce",srce_pt);

    multi1d<int> latsize;
    read(input_smeared_prop_xml,"/SinkSmear/Config_info/Params/HMCTrj/nrow",latsize);

    std::string traj;
    read(input_smeared_prop_xml,"/SinkSmear/Config_info/Params/MCControl/StartUpdateNum",traj);

    std::string srcesmear;
    read(input_smeared_prop_xml,"/SinkSmear/PropSource/Source/SourceType",srcesmear);

    std::string sinksmear;
    read(input_smeared_prop_xml,"/SinkSmear/PropSink/Sink/SinkType",sinksmear);

    std::string basename=""; //;+srcesmear+"_"+sinksmear+"_"+traj;
    std::string dirname=params.filedir;


    PASS(268);
    
    /////////////// New stuff ///////////////
    XMLReader xml_in;
      Cfg_t cfg;
    Real mass = 0.0;
    multi1d<int> bndry_cndtns;
    int ncg_had = 0;
    int j_decay = 0;
    multi1d<int> t_srce;
    multi1d<int> nrow;
    
    int dir;
    int steps;
    try
      {
	xml_in.open(getXMLInputFileName());
      }
    catch (...)
      {
	cerr << "Error in Reading input XML:" << endl;
	QDP_abort(1);
      }
    try
      {
	read(xml_in, "/chroma/Param/nrow", nrow);
	read(xml_in, "/chroma/Cfg", cfg);
	read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/FermionAction/Mass", mass);
	read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/FermionAction/FermState/FermionBC/boundary", bndry_cndtns);
	read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/j_decay", j_decay);
	read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/t_srce", t_srce);
	read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Displacement/dir", dir);
	read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Displacement/steps", steps);
      }
    catch (const string &e)
      {
	QDPIO::cerr << "Parsing XML: " << e << endl;
	QDP_abort(1);
      }

    PASS(309);
    
    Layout::setLattSize(nrow);
  Layout::create();

  PASS(314);

  LatticePropagator quark_prpgtr = zero;
  LatticePropagator quark_src = zero;
  quark_src = getPointSource(t_srce, j_decay);

  XMLReader rd(xml_in, "/chroma/Param/InlineMeasurements2/elem/propagator/Param");

  SimpleFermBCParams bc_prmtrs;

  GroupXML_t InvertParam = readXMLGroup(rd, "InvertParam", "invType");

  //boundary conditions
  bc_prmtrs.boundary.resize(Nd);
  for (int i = 0; i < Nd - 1; i++)
  {
    bc_prmtrs.boundary[i] = bndry_cndtns[i]; //BC_TYPE_PERIODIC;
  }
  bc_prmtrs.boundary[Nd - 1] = bndry_cndtns[Nd - 1]; //BC_TYPE_ANTIPERIODIC;
  

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

  //rgauge(u);

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
    push(xml_out, "meson_crrltr");

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
//drop gamma(8) for test: zero dispalcement FH_prpgtr should be identical to quark prpgtr
//    FH_src = gaugelink * ( quark_prpgtr_Delta );


    //======================Feynman-Hellman Source=========================

    Handle<FermState<T, P, P>> state(S.createState(u));

    LinOpSysSolverCGEnv::registerAll();

    QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;

    //===============2nd inversion Feynman-Hellman propagator==================

    XMLFileWriter xml_out("fh_xmldat");
    push(xml_out, "FH_meson_crrltr");
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
        ////////////////////////////////////////




    

    int Nt=latsize[3];

    // Boilerplate stuff to the output XML
    push(xml_out, "my_measurement");
    write(xml_out, "update_no", update_no);

    // Write info about the program
    proginfo(xml_out);

    // Write out the input
    params.write(xml_out, "Input");

    SpinMatrix cg5P_pstv = Cg5();
  
  //!!!! IMPORTNT
    SpinMatrix T_unpol = Tunpol();
    
    LatticeComplex d_tr_prpgtrs_1 =
      trace(T_unpol * traceColor(quark_prpgtr * traceSpin(quarkContract13(FH_quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr))));

    LatticeComplex d_tr_prpgtrs_2 =
      trace(T_unpol * traceColor(quark_prpgtr * quarkContract13(FH_quark_prpgtr * cg5P_pstv, cg5P_pstv * quark_prpgtr)));

    LatticeDouble d_tr_prpgtrs_Re = real(d_tr_prpgtrs_1 + d_tr_prpgtrs_2);
    LatticeDouble d_tr_prpgtrs_Im = imag(d_tr_prpgtrs_1 + d_tr_prpgtrs_2);


    /* First things first:  Need to convert prop from DeGrand-Rossi to Dirac-Pauli */
    const SpinMatrixD U =  PauliToDRMat();
    LatticePropagator prop_PDu = adj(U) * quark_prpgtr * U;
    LatticePropagator prop_PDd = adj(U) * FH_quark_prpgtr * U;
  
    SpinMatrix cg5P_pstv_PlDrc = adj(U) * Cg5() * U;
    SpinMatrix T_unpol_PlDrc = adj(U) * Tunpol() * U;


    LatticeComplex d_tr_prpgtrs_PlDrc_1 =
      trace(T_unpol_PlDrc * traceColor(prop_PDu * traceSpin(quarkContract13(prop_PDd * cg5P_pstv_PlDrc, cg5P_pstv_PlDrc * prop_PDu))));

    LatticeComplex d_tr_prpgtrs_PlDrc_2 =
      trace(T_unpol_PlDrc * traceColor(prop_PDu * quarkContract13(prop_PDd * cg5P_pstv_PlDrc, cg5P_pstv_PlDrc * prop_PDu)));

    LatticeDouble d_tr_prpgtrs_PlDrc_Re = real(d_tr_prpgtrs_PlDrc_1 + d_tr_prpgtrs_PlDrc_2);
    LatticeDouble d_tr_prpgtrs_PlDrc_Im = imag(d_tr_prpgtrs_PlDrc_1 + d_tr_prpgtrs_PlDrc_2);

    multi1d<int> pN(4);

//==================out put quark_prpgtr, FH_quark_prpgtr ==================//
LatticePropagatorReImPrt FH_quark_prpgtrRe = real(FH_quark_prpgtr);
LatticePropagatorReImPrt FH_quark_prpgtrIm = imag(FH_quark_prpgtr);
LatticePropagatorReImPrt quark_prpgtr_Re = real(quark_prpgtr);
LatticePropagatorReImPrt quark_prpgtr_Im = imag(quark_prpgtr);
HDF5Writer proton_PDF_file;
proton_PDF_file.open("proton_PDF_Tom.h5", HDF5Base::ate);
proton_PDF_file.set_stripesize(1048576);
  proton_PDF_file.push("FH_quark_prpgtr");
    proton_PDF_file.write("Re", FH_quark_prpgtrRe, HDF5Base::trunc);
    proton_PDF_file.write("Im", FH_quark_prpgtrIm, HDF5Base::trunc);
  proton_PDF_file.pop();
  proton_PDF_file.push("quark_prpgtrs");
    proton_PDF_file.write("Re", quark_prpgtr_Re, HDF5Base::trunc);
    proton_PDF_file.write("Im", quark_prpgtr_Im, HDF5Base::trunc);
  proton_PDF_file.pop();
  proton_PDF_file.push("d_tr_prpgtrs_Tom");
    proton_PDF_file.write("Re", d_tr_prpgtrs_Re, HDF5Base::trunc);
    proton_PDF_file.write("Im", d_tr_prpgtrs_Im, HDF5Base::trunc);
  proton_PDF_file.pop();
    proton_PDF_file.push("d_tr_prpgtrs_PlDrc_Tom");
    proton_PDF_file.write("Re", d_tr_prpgtrs_PlDrc_Re, HDF5Base::trunc);
    proton_PDF_file.write("Im", d_tr_prpgtrs_PlDrc_Im, HDF5Base::trunc);
  proton_PDF_file.pop(); 
proton_PDF_file.close();


    pN[0]=0;pN[1]=0;pN[2]=0;pN[3]=0;
    
    
    QDPIO::cout << "Propagator ID is " << params.smeared_prop_id << endl;

    /* instantiate all CalLAT classes here */
     
    Nucleon_system nucleon;
    /*
      Lambda_system lambda;
      Cascade_system cascade;
      Delta_system delta;
    */
    //    Omega_system omega;

    QDPIO::cout << params.nucleon_system <<" Beginning " << endl;
    nucleon.doChecks = params.consistencyTests;
    nucleon.initialize(srce_pt,pN,latsize,dirname,basename,"Avg");
    //    nucleon.calcDirect(prop_PD);
    //    nucleon.writeMass();
//    nucleon.calcExchange(prop_PD,prop_PD);
//for 3 point function
//    nucleon.calcBlock(prop_PDu,prop_PDd,112);
//for 2 point function
     nucleon.calcBlock(prop_PDu,prop_PDu,112);    

    nucleon.calcMassFromBlocks();
    //    nucleon.calc_G1pG1p_A1p();
    nucleon.free();


    pop(xml_out);
    measure_time.stop();
    QDPIO::cout << params.nucleon_system << ": total time = "
		<< measure_time.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << params.nucleon_system << ": ran successfully" << endl;
    END_CODE();
  } 

};

