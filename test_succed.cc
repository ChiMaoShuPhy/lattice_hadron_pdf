/*! \file
 * \brief Main code for tutorial 2
 */

#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"

using namespace std;
using namespace QDP;
using namespace Chroma;



//
// Main program
//
int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);
  
  // Put this in to enable profiling etc.
  START_CODE();

  // Instantiate xml reader for DATA
  // if you used -i input file, Chroma::getXMLInputFileName() returns
  //  the filename you specified
  // otherwise it will return "DATA"
  XMLReader xml_in;

  try { 
    xml_in.open(Chroma::getXMLInputFileName());
  }
  catch(...) { // Catch ALL Exceptions
    QDPIO::cerr << "Problem reading input XML file: " 
		<< Chroma::getXMLInputFileName() 
		<< endl;
    QDP_abort(1);  // Semi Graceful exit.
  }


  //1 D array to hold the size of the lattice 
  multi1d<int> nrow;
  //Chroma defined structure to hold the configuration format and file name
  Cfg_t cfg;

  //The source file name
  std::string src_file;  

  Real kappa = 0.;
  //Read the parameters from the input XML using absolute paths.
  try { 

    read(xml_in, "/tutorial2/Param/nrow", nrow);
    read(xml_in, "/tutorial2/Cfg", cfg);
    read(xml_in, "/tutorial2/Src/src_file", src_file);
    read(xml_in, "/tutorial2/Inversion/Kappa", kappa);

  }
  catch(const string& e) { 
    QDPIO::cerr << "Parsing XML: " << e << endl;
    QDP_abort(1);  // Semi Graceful exit.
  }


  // Terminal I/O to stdout.
  QDPIO::cout << "Tutorial 2" << endl;

  // Set lattice size, shape, etc.
  // using the nrow param from the Prameter structure
  Layout::setLattSize(nrow);

  // Initialise
  Layout::create();
  
  // -------------------------------------------------------------------
  // Read in the configuration along with relevant information.
  // -------------------------------------------------------------------

  // The gauge field may contain some file XML and other 
  // information (gauge XML): eg in the future the ILDG header(?)
  // We create empty XML readers for this, the gauge reading routine
  // will fill them in
  XMLReader gauge_file_xml, gauge_xml;

  // The gauge field itself
  multi1d<LatticeColorMatrix> u(Nd);
  
  // A convenience routine to read various gauge field formats
  gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);
  // ---------------- Gauge Reading Done -------------------------------


  // ---------------- Read in a source -----------------------------

  // Propagator also contains headers: a file XML and a record XML
  // as well as the binary data -- as dictated by SciDAC standards.
  // The files are expected to be in SciDAC format!
  XMLReader src_file_xml, src_record_xml;

  // The propagator itself
  LatticePropagator quark_source;
 
  QDPFileReader src_reader(src_file_xml,
			   src_file,     
			   QDPIO_SERIAL);

  src_reader.read(src_record_xml, quark_source);

  src_reader.close();

  // --------------- Source reader done ---------------------------


  // --------------- Insert some new code below -----------------------

  try{

    XMLReader rd(xml_in,"/tutorial2/Inversion");
    GroupXML_t  invParam =  readXMLGroup(rd, "InvertParam", "invType");

    typedef QDP::LatticeFermion T;
    typedef QDP::multi1d<LatticeColorMatrix> P;

    int ncg_had = 0;

    const Real mass = kappaToMass(kappa);
    QDPIO::cout << "kappa "<<kappa<<" mass "<<mass<<endl;

  // Parameters for boundary conditions
    SimpleFermBCParams bc_params;

    bc_params.boundary.resize(Nd);
    for (int i=0; i < Nd-1; i++){
      bc_params.boundary[i] = BC_TYPE_PERIODIC;
    }
    bc_params.boundary[Nd-1] = BC_TYPE_ANTIPERIODIC;
    
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    
  // The boundary conditions themselves - we will get one from the
  // heap using new and drop it into the handle, since that is what
  // the CreateFermState constructor wants.
    Handle<FermBC<T,P,P> > fbc_handle( new SimpleFermBC< T,P,P >(bc_params)) ;
    
  // Create a FermState Creator with boundaries
    Handle<CreateFermState<T,P,P> > cfs (new CreateSimpleFermState<T,P,P>(fbc_handle));
    
  // Now create the FermAct using the FermState Creator
    EvenOddPrecWilsonFermAct S(cfs, mass);
    
  // Use the FermAct to create a linear Operator over the FermState
    Handle< FermState<T,P,P> > state( S.createState(u) );
    
    LatticePropagator quark_propagator = zero;
    
    XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

    push(xml_out, "tutorial2");
    
    LinOpSysSolverCGEnv::registerAll();
    
    QuarkSpinType quarkSpinType = QUARK_SPIN_TYPE_FULL;
    
    S.quarkProp(quark_propagator,
		xml_out,
		quark_source,
		state,
		invParam,
		quarkSpinType,
		ncg_had);
    
    xml_out.flush();
    
    pop(xml_out);
  }
  catch (const std::string &e){
    QDPIO::cout << "Error trying to perform the inversion: "<< e<<endl;
    QDP_abort(1);
  }
  
  
  END_CODE(); // End 


  Chroma::finalize();

  exit(0);
}
