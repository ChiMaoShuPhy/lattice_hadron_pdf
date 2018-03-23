/* 
 * File:   nucleon_meas.h
 * Author: tomluu
 *
 * Created on December 6, 2013, 10:07 AM
 */

#ifndef NUCLEON_MEAS_H
#define	NUCLEON_MEAS_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "nucleon.h"

namespace Chroma 
{ 
  // A namespace for this particular   measurement
  namespace NucleonEnv 
  {
    extern const std::string name;
    bool registerAll();
  

    //! Parameter structure
    struct Params 
    {
      // Default constructor -- forward declaration
      Params();

      // Construct from XML -- forward declaration
      Params(XMLReader& xml_in, const std::string& path);

      // Write myself out
      void write(XMLWriter& xml_out, const std::string& path);

      // How often should I measure me in an HMC run
      unsigned long frequency;

      // Various parameters taken from XML file
      std::string nucleon_system;
      std::string gauge_id;
      std::string smeared_prop_id;
      std::string filedir;
      bool consistencyTests; /* not mandatory */

      std::string xml_file; // Support output to own private XML File
    }; // struct
  }; // namespace NucleonEnv

  class Nucleon : public AbsInlineMeasurement 
  {
  public:
    // Constructor -- default -- do nothing
    ~Nucleon() {}

    // Constructor -- from param struct -- copy param struct inside me
    Nucleon(const NucleonEnv::Params& p) : params(p) {}

    // Constructor -- copy constructor -- copy param struct from argument
    Nucleon(const Nucleon& p) : params(p.params) {}

    // Boiler plate
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 


  private:

    void func(const unsigned long update_no,XMLWriter& xml_out);
    
    NucleonEnv::Params params;
  };

};

#endif	/* NUCLEON_MEAS_H */

