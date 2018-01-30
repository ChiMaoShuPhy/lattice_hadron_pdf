#ifndef GAUGE_LINK_H
#define GAUGE_LINK_H

LatticeColorMatrix gauge_shift(const multi1d<LatticeColorMatrix> &U, const LatticeColorMatrix &v, const int sign, const int dir)
{
    if (sign == 0)
        return v;
    return sign > 0 ? shift(U[dir] * v, +1, dir)
                    : shift(adj(U[dir]) * v, -1, dir);
}
// create a gauge link in \bar psi(n)L(n,n+\Delta\hat \mu)psi(n+\Delta \hat \mu)
// L(n,n+\Delta\hat\mu)=U_\mu(n)*U_\mu(n+\hat\mu)*...*U_\mu(n+(\Delta-1)\hat\mu)
LatticeColorMatrix gauge_link(const multi1d<LatticeColorMatrix> &U, const int sign, const int dir, const int displacement)
{
    Real one = 1.0;

    LatticeColorMatrix gaugelink = one;

    if (displacement == 0)
    {
        return one;
    }
    else if (displacement > 0)
    {
        for (int stp = 1; stp < displacement; stp++)
        {
            gaugelink = gauge_shift(U, gaugelink, 1, dir);
        }
        return U[dir] * gaugelink;
        //  return gaugelink;
    }
    else
    {
        for (int stp = -1; stp >= displacement; stp--)
        {
            gaugelink = gauge_shift(U, gaugelink, -1, dir);
        }
        return gaugelink;
    }
}


#endif