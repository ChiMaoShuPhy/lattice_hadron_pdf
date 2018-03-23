#include "kn_meas.h"
#include "chroma.h"
#include "spin_matrix.h"
#include <iostream>
#include <string>
#include "propstuff.h"

using namespace QDP;

namespace Chroma 
{ 
  // Fill in the blanks
  namespace KaonNucleonEnv 
  { 

    // Function to register with a factory
    // This is boilerplate stuff
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      // Create a Params from the XML
      // Use it to create a KaonNucleon
      return new KaonNucleon( Params(xml_in, path) );
    }

    // The name of my measurement for the XML file
    // Change this for each measurement
    const std::string name = "KN_SYSTEM";

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
    Params::Params() { kn_system="";frequency = 0; gauge_id=""; smeared_nucleon_q=""; xml_file=""; filedir=""; consistencyTests = false;}
    
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
	  read(paramtop, "Name", kn_system);

	  // Read in the file directory
	  read(paramtop,"Param/file_dir",filedir);
          
	  // Read in the output propagator/source configuration info
	  read(paramtop, "NamedObject/smeared_nucleon_q", smeared_nucleon_q);

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
      QDP::write(xml_out, "smeared_nucleon_q", smeared_nucleon_q);

      QDP::write(xml_out, "gauge_id", gauge_id);

      if( xml_file != "" ){ 
	QDP::write(xml_out, "xml_file", xml_file);
      }


      pop(xml_out);
    }
  }

  void 
  KaonNucleon::operator()(long unsigned int update_no,
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

  void 
  KaonNucleon::func(unsigned long int update_no, XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch measure_time;
    measure_time.reset();
    measure_time.start();


    // Test that the gauge configuration and the propagator we need
    // exist in the map.
    XMLBufferWriter gauge_xml;
    XMLBufferWriter smeared_prop_xml;
    try
      {
	// Try and get at the gauge field if it doesn't exist 
	// an exception will be thrown.
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);
	TheNamedObjMap::Instance().get(params.gauge_id).getRecordXML(gauge_xml);

	// Do the same with the propagator 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.smeared_nucleon_q);

	// Get its record XML
	TheNamedObjMap::Instance().get(params.smeared_nucleon_q).getRecordXML(smeared_prop_xml);

      }
    catch( std::bad_cast ) 
      {

	// If the object exists and is of the wrong type we will 
	// land in this catch.
	QDPIO::cerr << params.kn_system << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
    catch (const string& e) 
      {
	// If the object is not in the map we will land in this 
	// catch
	QDPIO::cerr << params.kn_system << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }

    // If we got here, then the gauge field and prop are both in 
    // the map. Their XML will have been captured.
    // Let us bind the references to a local name so 
    // we don't have to type the long lookup string every time.
    //
    // Note const here means we can't change the field
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);


    const LatticePropagator& smeared_prop = 
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.smeared_nucleon_q);

    XMLReader input_smeared_prop_xml;
    TheNamedObjMap::Instance().get(params.smeared_nucleon_q).getRecordXML(input_smeared_prop_xml);

    multi1d<int> srce_pt;
    read(input_smeared_prop_xml,"/SinkSmear/PropSource/Source/t_srce",srce_pt);

    multi1d<int> latsize;
//    read(input_smeared_prop_xml,"/SinkSmear/Config_info/Params/HMCTrj/nrow",latsize);
    latsize.resize(4); //// Modified for Thorsten's code here!!!!!!!!!!!
    latsize[0]=4;latsize[1]=4;latsize[2]=4;latsize[3]=4;  //// Modified for Thorsten's code here!!!!!!!!!!!

    std::string traj;
//    read(input_smeared_prop_xml,"/SinkSmear/Config_info/Params/MCControl/StartUpdateNum",traj);
    traj = "separate"; //// Modified for Thorsten's code here!!!!!!!!!!!
    
    std::string srcesmear;
    read(input_smeared_prop_xml,"/SinkSmear/PropSource/Source/SourceType",srcesmear);
    srcesmear = "";
    
    std::string sinksmear;
    read(input_smeared_prop_xml,"/SinkSmear/PropSink/Sink/SinkType",sinksmear);
    sinksmear = "";
    
    std::string basename="_"+srcesmear+"_"+sinksmear+"_"+traj;
    std::string dirname=params.filedir;
    int Nt=latsize[3];

    // Boilerplate stuff to the output XML
    push(xml_out, "my_measurement");
    write(xml_out, "update_no", update_no);

    // Write info about the program
    proginfo(xml_out);

    // Write out the input
    params.write(xml_out, "Input");

    /* First things first:  Need to convert prop from DeGrand-Rossi to Dirac-Pauli */
    const SpinMatrixD U =  PauliToDRMat();

    /*
     Here is temporary code to check with Thorsten
     */
    LatticePropagator nucleonQ,kaonQ,kaonS;    
    loadProp("/Users/tomluu/Research/lqcd/callat/kn/props/prop.lat.tst.plain.0.u_p.hi.dat", nucleonQ);
    loadProp("/Users/tomluu/Research/lqcd/callat/kn/props/prop.lat.tst.plain.0.u_k.hi.dat", kaonQ);
    loadProp("/Users/tomluu/Research/lqcd/callat/kn/props/prop.lat.tst.plain.0.s_k.hi.dat", kaonS);
    nucleonQ = transpose(U) * nucleonQ * U;
    kaonQ = transpose(U) * kaonQ * U;
    kaonS = transpose(U) * kaonS * U;
    
//    LatticePropagator prop_PD;
//    prop_PD = transpose(U) * smeared_prop * U;

//    QDPIO::cout << "Propagator ID is " << params.smeared_nucleon_q << endl;

    /* instantiate kaon_nucleon class here */
     
    Kaon_Nucleon_system kn;
    multi1d<int> mom(4);

    mom[0]=mom[1]=mom[2]=mom[3]=0;
    
    
    
    QDPIO::cout << params.kn_system <<" Beginning " << endl;
    kn.doChecks = params.consistencyTests;
    kn.initialize(srce_pt,mom,srce_pt,mom,latsize,dirname,basename);
    kn.calcDirect(nucleonQ,kaonQ,kaonS);
    kn.calcExchange(nucleonQ,kaonQ,kaonS);
//    kn.calcDirect(prop_PD, prop_PD, prop_PD);
//    kn.calcExchange(prop_PD, prop_PD, prop_PD);
    kn.writeCorrelators();
    kn.free();

    pop(xml_out);
    measure_time.stop();
    QDPIO::cout << params.kn_system << ": total time = "
		<< measure_time.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << params.kn_system << ": ran successfully" << endl;
    END_CODE();
  } 

};
