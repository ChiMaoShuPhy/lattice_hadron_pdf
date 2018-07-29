#ifndef TIMESLICE_FT_H/*define*/
#define TIMESLICE_FT_H/*define*/

#include "chroma.h"

#define pi M_PI

using namespace std;
using namespace QDP;
using namespace Chroma;


//struct multi1d_int {multi1d<int>int4(4);};


template <typename Typ_rtn, typename Typ_vrbl>
Typ_rtn timesilce_FT(Typ_vrbl lattice_data,multi1d<int> &nrow_p, multi1d<int> &momentum_p, int j_decay, int Nd)
{
  LatticeDouble r_dot_p = zero;
  LatticeDouble coord_in_dir;

  multi1d<double> momentum(4);
for(int dir = 0; dir <Nd; dir++)
{ 
    momentum[dir] = momentum_p[dir]*2*pi/nrow_p[dir]; 
}
  for(int dir_i = 0; dir_i<Nd-1; dir_i++)
    {
     coord_in_dir = (LatticeDouble) Layout::latticeCoordinate(dir_i);
     r_dot_p = r_dot_p + momentum[dir_i]*coord_in_dir;  
    }

  LatticeComplex lattice_data_FT = (cos(r_dot_p) 
                         - cmplx(Real(0.0),Real(1.0))*sin(r_dot_p))*lattice_data;

  Set timeslice;
  timeslice.make(TimeSliceFunc(j_decay));
  int length = timeslice.numSubsets();
  Typ_rtn data_FT(length);
  data_FT = sumMulti(lattice_data_FT, timeslice);

  return data_FT;
}



#endif ///*define*/
