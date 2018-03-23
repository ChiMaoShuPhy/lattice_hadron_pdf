#include "permutations.h"
using namespace std ;

string sign2string(const Permutation& p){
  string s ;
  if(p.sign()>0)
    s = "+";
  else
    s = "-";
  return s ;
}


int epsilon(const multi1d<int>& sigma){
  int num(1),dnum(1) ;
  for(int i(0);i<sigma.size();i++)
    for(int j(i+1);j<sigma.size();j++)
      {
	num *= (sigma[j] - sigma[i]);
	dnum *= (j - i) ;
      }
  return num/dnum ;
}

int epsilon(int a, int b, int c){
  multi1d<int> i(3);
  i[0]=a ; i[1]=b; i[2]=c ;
  return epsilon(i);
}

int delta(int a, int b){
  if(a==b) return 1 ;
  return 0 ;
}
