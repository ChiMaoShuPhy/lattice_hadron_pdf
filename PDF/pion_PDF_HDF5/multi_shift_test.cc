#include <string>
#include "chroma.h"
#include "multi_shift.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

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
        //  W[z+Delta_z,z], Delta_z: steps
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

    LatticeInt coord_in_dir = Layout::latticeCoordinate(dir);
    LatticeInt shftd_coord_in_dir = multi_shift<LatticeInt>(coord_in_dir, steps, dir, nrow[dir], +1);

    //      LatticeDouble shft_phs_in_dir = LatticeRealD(pow(-1, (shftd_coord_in_dir + steps) / N));
    LatticeDouble shft_phs_in_dir = one;

    if(steps > 0)
    {
    shft_phs_in_dir = LatticeRealD(pow(-1,LatticeInt(shftd_coord_in_dir < steps)));
    }
    else if (steps < 0)
    {
    shft_phs_in_dir = LatticeRealD(pow(-1,LatticeInt(shftd_coord_in_dir >= nrow[dir] + steps)));
    }
    else
    {
      ;   
    }

    HDF5Writer multi_shift_file;
    multi_shift_file.open("multi_shift.h5", HDF5Base::ate);
    multi_shift_file.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB.

    multi_shift_file.write("coord_in_dir", LatticeRealD(coord_in_dir), HDF5Base::trunc);
    multi_shift_file.write("shftd_coord_in_dir", LatticeRealD(shftd_coord_in_dir), HDF5Base::trunc);
    multi_shift_file.write("shft_phs_in_dir", shft_phs_in_dir, HDF5Base::trunc);


    multi_shift_file.close();

    END_CODE();
    finalize();
    exit(0);
}
