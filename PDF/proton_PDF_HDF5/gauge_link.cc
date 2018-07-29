#include <string>
#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

Real one = 1.0;
Real zero = 0.0;

LatticeColorMatrix gauge_shift(const multi1d<LatticeColorMatrix> &U,
                               const LatticeColorMatrix &v,
                               const int sign, const int dir)
{
    if (sign == 0)
        return v;
    return sign > 0 ? shift(adj(U[dir]) * v, -1, dir)
                    : U[dir] * shift(v, +1, dir);
}

LatticeColorMatrix gauge_link(const multi1d<LatticeColorMatrix> &U, const int dir, const int displacement)
{
    LatticeColorMatrix gaugelink = one;

    if (displacement == 0)
    {
        gaugelink = one;
    }

    else if (displacement > 0)
    {
        for (int stp = 1; stp <= displacement; stp++)
        {
            gaugelink = gauge_shift(U, gaugelink, 1, dir);
        }
    }

    else
    {
        for (int stp = -1; stp >= displacement; stp--)
        {
            gaugelink = gauge_shift(U, gaugelink, -1, dir);
        }
    }
    return gaugelink;
}

int main(int argc, char **argv)
{
    initialize(&argc, &argv);
    START_CODE();


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
    int ncg_had = 0;
    int j_decay = 0;
    multi1d<int> t_srce;

    int dir;
    int steps;

    //read in
    try
    {
        read(xml_in, "/chroma/Param/nrow", nrow);
        //    cout<<"read nrow done"<<endl;
        read(xml_in, "/chroma/Cfg", cfg);
        //    cout<<"read cfg done"<<endl;
        //    cout<<"read mass done"<<endl;
        //    cout<<"read boundry condition done"<<endl;
        read(xml_in, "/chroma/Param/InlineMeasurements/elem/Param/Source/Source_Smearing/j_decay", j_decay);
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

    //boundary conditions
    XMLReader gauge_file_xml, gauge_xml;
    multi1d<LatticeColorMatrix> u(Nd);

    //HotSt(u);

    gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);

    unitarityCheck(u);

    LatticeColorMatrix gauge_link_orgnl = one;
    LatticeColorMatrix gauge_link_trnsfrmd = one;
    LatticeColorMatrix gauge_link_dffrnc = zero;

    LatticeColorMatrix g1 = one;
    LatticeColorMatrix g2 = one;
    LatticeColorMatrix g0 = one;

    gauge_link_orgnl = gauge_link(u, dir, steps);

    rgauge(u, g1);

    g0 = g1;

    if (steps == 0)
    {
        g2 = shift(g0, -1, dir);
    }
    else if (steps > 0)
    {
        for (int i = 1; i <= steps; i++)
        {
            g2 = shift(g0, -1, dir);
            g0 = g2;
        }
    }
    else
    {
        for (int i = -1; i >= steps; i--)
        {
            g2 = shift(g0, +1, dir);
            g0 = g2;
        }
    }

    gauge_link_trnsfrmd = gauge_link(u, +1, dir, steps);

    gauge_link_dffrnc = adj(g1) * gauge_link_trnsfrmd * g2 -
                        gauge_link_orgnl;

    XMLFileWriter xml_out("gauge_link_difference");
    push(xml_out, "difference");

    write(xml_out, "delta", gauge_link_dffrnc);

    xml_out.flush();
    pop(xml_out);

    END_CODE();
    finalize();
    exit(0);
}