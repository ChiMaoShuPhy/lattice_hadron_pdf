#ifndef __KaonNucleon_h__
#define __KaonNucleon_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "kn_sys.h"

namespace Chroma 
{ 
  // A namespace for this particular   measurement
  namespace KaonNucleonEnv 
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
      std::string kn_system;
      std::string gauge_id;
      std::string smeared_nucleon_q;
      std::string filedir;
      bool consistencyTests; /* not mandatory */

      std::string xml_file; // Support output to own private XML File
    }; // struct
  }; // namespace KaonNucleonEnv

  class KaonNucleon : public AbsInlineMeasurement 
  {
  public:
    // Constructor -- default -- do nothing
    ~KaonNucleon() {}

    // Constructor -- from param struct -- copy param struct inside me
    KaonNucleon(const KaonNucleonEnv::Params& p) : params(p) {}

    // Constructor -- copy constructor -- copy param struct from argument
    KaonNucleon(const KaonNucleon& p) : params(p.params) {}

    // Boiler plate
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 


  private:

    void func(const unsigned long update_no,XMLWriter& xml_out);
    
    KaonNucleonEnv::Params params;
  };

};

#endif
