#ifndef __multiD_h__
#define __multiD_h__

#include "qdp.h"
using namespace std;
using namespace QDP;


template<class T> class multi6d
{
 public:
  // Basic constructor
  multi6d() {F=0;n1=n2=n3=n4=n5=n6=sz=0;copymem=false;}
  multi6d(T *f,int ns6, int ns5, int ns4, int ns3, int ns2, int ns1) {F=f; n1=ns1; n2=ns2; n3=ns3; n4=ns4; n5=ns5; n6=ns6; sz=n1*n2*n3*n4*n5*n6; copymem=true;}
  explicit multi6d(int ns6, int ns5, int ns4, int ns3, int ns2, int ns1) {copymem=false;F=0;resize(ns6,ns5,ns4,ns3,ns2,ns1);}
  ~multi6d() {if (! copymem) {delete[] F;}}
  
  //! Copy constructor
  multi6d(const multi6d& s): copymem(false), n1(s.n1), n2(s.n2), n3(s.n3), n4(s.n4), n5(s.n5),n6(s.n6), sz(s.sz), F(0)
    {
      resize(n6,n5,n4,n3,n2,n1);
      
      for(int i=0; i < sz; ++i)
	F[i] = s.F[i];
    }
    
    //! Allocate mem for the array 
  void resize(int ns6, int ns5, int ns4, int ns3, int ns2, int ns1) 
    {
      if(copymem) {
	QDPIO::cerr <<"multi6d block: invalid resize of a copy of memory" << endl;
	exit(1);
      }
      
      // Only delete if the array is not NULL. If it is NULL
      // deleting may be bad
      if ( F != 0x0 ) {
	delete[] F; 
      }
      
      n1=ns1; n2=ns2; n3=ns3; n4=ns4; n5=ns5; n6=ns6; sz=n1*n2*n3*n4*n5*n6; F = new(nothrow) T[sz];
      if( F == 0x0 ) { 
	QDP_error_exit("Unable to new memory in multi6d block::resize(%d,%d,%d,%d,%d,%d)\n",ns6,ns5,ns4,ns3,ns2,ns1);
      }
    }
    
  //! Size of array
  int size1() const {return n1;}
  int size2() const {return n2;}
  int size3() const {return n3;}
  int size4() const {return n4;}
  int size5() const {return n5;}
  int size6() const {return n6;}

  // Return reference to an element
  T& operator()(int n, int m, int l, int k, int j, int i) {return F[i+n1*(j+n2*(k+n3*(l+n4*(m+n5*(n)))))];}
  const T& operator()(int n, int m, int l, int k, int j, int i) const {return F[i+n1*(j+n2*(k+n3*(l+n4*(m+n5*(n)))))];}

  // Return multi4d reference
  multi5d<T> operator[](int n) {return multi5d<T>(F+n1*n2*n3*n4*n5*n,n5,n4,n3,n2,n1);}
  const multi5d<T> operator[](int n) const {return multi5d<T>(F+n1*n2*n3*n4*n5*n,n5,n4,n3,n2,n1);}

private:
  bool copymem;
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int sz;
  T *F;
};


#endif
