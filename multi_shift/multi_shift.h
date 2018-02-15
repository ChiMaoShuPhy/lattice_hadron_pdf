/*==========================================
Multiple shift with boundary condition imposed.
multi_shit<typ>(typ objct, int steps , int dir, int N, int bndry_cndtn)
returns the shifted objct.

steps      : number of steps in unit of lattice spacing.
dir        : direction index
N          : number of total lattice spacing in direction dir
bndry_cndtn: +1: periodic boundary condition f(x+N) = f(x)
             -1: anti-periodic boundary condition f(x+N) = -f(x)
==========================================*/
#ifndef MULTI_SHIFT
#define MULTI_SHIFT

#include "chroma.h"

using namespace std;


template <typename Typ>
Typ multi_shift(Typ objct, int steps, int dir, int N, int bndry_cndtn)
{
    Typ objct_tmp1 = objct;
    Typ objct_tmp2 = objct;
    Real one = 1.0;
    Real zero = 0.0;
    if (bndry_cndtn == 1)//periodic boundary condition
    {
        if (steps > 0)
        {
            for (int stp = 1; stp <= steps; stp++)
            {
                objct_tmp2 = shift(objct_tmp1, +1, dir);
                objct_tmp1 = objct_tmp2;
            }
        }

        else if (steps < 0)
        {
            for (int stp = -1; stp >= steps; stp--)
            {
                objct_tmp2 = shift(objct_tmp1, -1, dir);
                objct_tmp1 = objct_tmp2;
            }
        }

        else
        {
            objct_tmp1 = objct;
        }

        return objct_tmp1;
    }
    else if (bndry_cndtn == -1)//anti-periodic boundary condition
    {
        LatticeInt coord_in_dir = Layout::latticeCoordinate(dir);//cretaed indcies of each lattice site in diraction dir;
        LatticeInt shftd_coord_in_dir = multi_shift<LatticeInt>(coord_in_dir, steps, dir, N, +1);
        //shifting indcies (steps) steps in diraction dir;

        LatticeDouble shft_phs_in_dir = one;//shft_phs_in_dir: the phase factor on each site

        if (steps > 0)
        {
            shft_phs_in_dir = LatticeRealD(pow(-1, LatticeInt(shftd_coord_in_dir < steps))); //those sites acrossing the bondary in direction dir, receives a -1 factor
            //e.g before shifting               W = 0  1  2  3...27 28 29 30 31
            //   after shifting(steps=3,dir =3) W'= 3  4  5  6...30 31 0  1  2
            //QDP shift: W'=shift(W,+1,dir) gives W'[n] = W[n+1]
            //W'[29;;31] = 0 1 2, have accrossed boundary dues receives a -1 factor,  (-1)^(W'[n]<steps) 
            // where N=32 in the example, W dennotes coord_in_dir and W' denotes shftd_coord_in_dir in the code
        }
        else if (steps < 0)
        {
            shft_phs_in_dir = LatticeRealD(pow(-1, LatticeInt(shftd_coord_in_dir >= N + steps)));
            //e.g before shifting                W = 0  1  2  3...27 28 29 30 31
            //   after shifting(steps= -3,dir=3) W'= 29 30 31 0...30 31 0  27 28
            //QDP shift: W'=shift(W,-1,dir) gives W'[n] = W[n-1]
            //W'[0;;3] = 29 30 31 have accrossed boundary dues receives a -1 factor,  (-1)^(W'[n]>= N + steps) 
            // where N=32 in the example, W dennotes coord_in_dir and W' denotes shftd_coord_in_dir in the code
        }
        else
        {
            ;
        }

        return shft_phs_in_dir * (multi_shift<Typ>(objct, steps, dir, N, +1));
    }

    else
    {
        QDPIO::cerr << "Invalid boundary condition for multi_shift " << endl;
        QDP_abort(1);
    }
}

#endif