
#include "chromabase.h"


#ifndef _permutations_h_
#define _permutations_h_

using namespace std;
using namespace QDP;


class Permutation {
private:
  multi1d<int> sigma ;

  int count ;
  int NumPerm ;

  bool check(const multi1d<int>& p){
   
    for(int i(0);i<p.size();i++)
      for(int j(i+1);j<p.size();j++)
	if(p[i]==p[j]) 
	  return false ;
    
    return true ;
  }

  bool bound_checks(const multi1d<int>& p){
    for(int i(0);i<p.size();i++)
      if(p[i]>=p.size())
	return false ;
    return true ;
  }

  void calc_NumPerm(){
    NumPerm = 1 ;
    for(int i(0);i<sigma.size();i++){
      NumPerm *= (i+1);
    }
  }

public:
  Permutation(int n){
    sigma.resize(n);
    for(int i(0);i<n;i++)
      sigma[i] = i ;
  }

  Permutation(const multi1d<int>& p){
    sigma.resize(p.size());
    set(p);
  }

  void set(const multi1d<int>& p){
    if(p.size() != sigma.size()){
      QDPIO::cerr<<"incompatible sizes. Abord!\n";
      QDP_abort(123);
    }
    if(bound_checks(p)){
      if(check(p)) 
	sigma = p ;
      else{
	QDPIO::cerr<<"not a permutation!\n";
	QDP_abort(112) ;
      }
    }
    else{
      QDPIO::cerr<<"not a permutation. out of bounds!\n";
      QDP_abort(124) ;
    }
  }
  
  int sign() const {
    int num(1),dnum(1) ;
    for(int i(0);i<sigma.size();i++)
      for(int j(i+1);j<sigma.size();j++)
	{
	  //QDPIO::cerr<<"( "<<sigma[j]<<" - "<<sigma[i]<<" )/( "<<j<<" - "<<i<<" )\n";
	  num *= (sigma[j] - sigma[i]);
	  dnum *= (j - i) ;
	}
    //QDPIO::cerr<<num/dnum<<endl ;
    return num/dnum ;
  }

  void begin(){
    count = 0 ;
    for(int i(0);i<sigma.size();i++)
      sigma[i] = i ;
    calc_NumPerm();
  }

  int operator[](int i){
    return sigma[i] ;
  }

  void next(){
    int n(sigma.size());
    
    int i(0) ;
  
    sigma[i]++ ;
    while(sigma[i]>=n)
      {
	sigma[i]=0 ;
      if(i<n-1)
        {
          i++ ;
          sigma[i]++ ;
        }
      }
    if(check(sigma)){
      count++ ;
      return ;
    }
    else
      next();
  }

  bool not_end(){
    return (count<NumPerm) ;
  }

  bool the_end(){
    return (count==NumPerm-1) ;
  }

  
};

string sign2string(const Permutation& p);
int epsilon(const multi1d<int>& sigma);
int epsilon(int a, int b, int c);
int delta(int a, int b);

#endif
