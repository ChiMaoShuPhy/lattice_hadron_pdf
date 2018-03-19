#ifndef __udd_h__
#define __udd_h__

#include "qdp.h"
#include "myMultiD.h"

using namespace std;
using namespace QDP;

template<class T> class udd
{
 public:
  // Basic constructor
  udd() {F=0;n1=n2=n3=n4=n5=n6=n7=sz=0;copymem=false;}
  udd(T *f, int ns7, int ns6, int ns5, int ns4, int ns3, int ns2, int ns1) {F=f; n1=ns1; n2=ns2; n3=ns3; n4=ns4; n5=ns5; n6=ns6; n7=ns7;  sz=n1*n2*n3*n4*n5*n6*n7; copymem=true;}
  explicit udd(int ns7, int ns6, int ns5, int ns4, int ns3, int ns2, int ns1) {copymem=false;F=0;resize(ns7,ns6,ns5,ns4,ns3,ns2,ns1);}
  ~udd() {if (! copymem) {delete[] F;}}
  
  //! Copy constructor
 udd(const udd& s): copymem(false), n1(s.n1), n2(s.n2), n3(s.n3), n4(s.n4), n5(s.n5), n6(s.n6), n7(s.n7), sz(s.sz), F(0)
    {
      resize(n7,n6,n5,n4,n3,n2,n1);
      
      for(int i=0; i < sz; ++i)
	F[i] = s.F[i];
    }
    
    //! Allocate mem for the array 
    void resize(int ns7, int ns6, int ns5,int ns4, int ns3, int ns2, int ns1) 
    {
      if(copymem) {
	cerr <<"udd block: invalid resize of a copy of memory" << endl;
	exit(1);
      }
      
      // Only delete if the array is not NULL. If it is NULL
      // deleting may be bad
      if ( F != 0x0 ) {
	delete[] F; 
      }
      
      n1=ns1; n2=ns2; n3=ns3; n4=ns4; n5=ns5; n6=ns6; n7=ns7; sz=n1*n2*n3*n4*n5*n6*n7; F = new(nothrow) T[sz];
      if( F == 0x0 ) { 
	QDP_error_exit("Unable to new memory in udd block::resize(%d %d %d %d,%d,%d,%d)\n",ns7,ns6,ns5,ns4,ns3,ns2,ns1);
      }
    }
    
  //! Size of array
  int size1() const {return n1;}
  int size2() const {return n2;}
  int size3() const {return n3;}
  int size4() const {return n4;}
  int size5() const {return n5;}
  int size6() const {return n6;}
  int size7() const {return n7;}

  // Return reference to an element
  T& operator()(int o, int n, int m, int l, int k, int j, int i) {return F[i+n1*(j+n2*(k+n3*(l+n4*(m+n5*(n+n6*(o))))))];}
  const T& operator()(int o, int n, int m, int l, int k, int j, int i) const {return F[i+n1*(j+n2*(k+n3*(l+n4*(m+n5*(n+n6*(o))))))];}
  //  const T& operator()(int l, int k, int j, int i) const {return F[i+n1*(j+n2*(k+n3*(l+n4*(m+n5*(n+n6*(o))))))];}


  // Return multi6d reference
  multi6d<T> operator[](int o) {return multi6d<T>(F+n1*n2*n3*n4*n5*n6*o,n6,n5,n4,n3,n2,n1);}
  const multi6d<T> operator[](int o) const {return multi6d<T>(F+n1*n2*n3*n4*n5*n6*o,n6,n5,n4,n3,n2,n1);}

  // Return multi1d reference
  //  multi1d<T> operator[][][][][][](int o, int n, int m, int l, int k, int j) {return multi1d<T>(F+n1*(j+n2*(k+n3*(l+n4*(m+n5*(n+n6*(o)))))),n1);}
  //  const multi1d<T> operator[][][][][][](int o, int n, int m, int l, int k, int j) const {return multi1d<T>(F+n1*(j+n2*(k+n3*(l+n4*(m+n5*(n+n6*(o)))))),n1);}

  //! Equal operator uses underlying = of T
  udd<T>& operator=(const udd<T>& s1)
    {
      resize(s1.size7(), s1.size6(), s1.size5(), s1.size4(), s1.size3(), s1.size2(), s1.size1());
      
      for(int i=0; i < sz; ++i)
	F[i] = s1.F[i];
      return *this;
    }

  template<class T1> udd<T>& operator=(const T1& s1)
    { /* sets equal to a constant */
      if(F==0)
	{
	  cerr<<"udd block: left hand side not initialized" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] = s1;
      return *this;
    }

  //! Add-replace operator uses underlying += of T
  udd<T>& operator+=(const udd<T>& s1)
    {
      if ( size1() != s1.size1() || size2() != s1.size2()
	   || size3() != s1.size3() || size4() != s1.size4()
	   || size5() != s1.size5() || size6() != s1.size6()
	   || size7() != s1.size7() )
	{
	  cerr << "udd block: Sizes incompatible in +=" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] += s1.F[i];
      return *this;
    }

  template<class T1> udd<T>& operator+=(const T1& s1)
    { /* adds a constant */
      if(F==0)
	{
	  cerr<<"udd block: left hand side not initialized" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] += s1;
      return *this;
    }

  //! Subtract-replace operator uses underlying -= of T
  udd<T>& operator-=(const udd<T>& s1)
    {
      if ( size1() != s1.size1() || size2() != s1.size2()
	   || size3() != s1.size3() || size4() != s1.size4()
	   || size5() != s1.size5() || size6() != s1.size6()
	   || size7() != s1.size7() )
	{
	  cerr << "udd block: Sizes incompatible in -=" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] -= s1.F[i];
      return *this;
    }

  template<class T1> udd<T>& operator-=(const T1& s1)
    { /* subtracts a constant */
      if(F==0)
	{
	  cerr<<"udd block: left hand side not initialized" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] -= s1;
      return *this;
    }

  //! Multiply-replace operator uses underlying *= of T
  udd<T>& operator*=(const udd<T>& s1)
    {
      if ( size1() != s1.size1() || size2() != s1.size2()
	   || size3() != s1.size3() || size4() != s1.size4()
	   || size5() != s1.size5() || size6() != s1.size6()
	   || size7() != s1.size7() )
	{
	  cerr << "udd block: Sizes incompatible in *=" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] *= s1.F[i];
      return *this;
    }

  template<class T1> udd<T>& operator*=(const T1& s1)
    { /* multiplies a constant */
      if(F==0)
	{
	  cerr<<"udd block: left hand side not initialized" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] *= s1;
      return *this;
    }
  //! Divide-replace operator uses underlying /= of T
  udd<T>& operator/=(const udd<T>& s1)
    {
      if ( size1() != s1.size1() || size2() != s1.size2()
	   || size3() != s1.size3() || size4() != s1.size4() 
	   || size5() != s1.size5() || size6() != s1.size6()
	   || size7() != s1.size7() )
	{
	  cerr << "udd block: Sizes incompatible in /=" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] /= s1.F[i];
      return *this;
    }

  template<class T1> udd<T>& operator/=(const T1& s1)
    { /* divides a constant */
      if(F==0)
	{
	  cerr<<"udd block: left hand side not initialized" << endl;
	  exit(1);
	};
      
      for(int i=0; i < sz; ++i)
	F[i] /= s1;
      return *this;
    }

  
  multi1d<Complex> uddTrace(int l, int k, int j) const
  {
    multi1d<Complex> traced(n1);
    
    for( int i=0; i < n1; ++i)
      traced[i] = trace( F[ i + n1*( j + n2*( k + n3*( l ) ) ) ] );

    return traced;
  }
  

private:
  bool copymem;
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int n7;
  int sz;
  T *F;
};

// add two udd blocks
template< typename T > inline udd<T> operator+(const udd<T>& a, const udd<T>& b)
{
  udd<T> c(a);
  c+=b;
  return c;
}

// subtract two udd blocks
template< typename T > inline udd<T> operator-(const udd<T>& a, const udd<T>& b)
{
  udd<T> c(a);
  c-=b;
  return c;
}

// multiply two udd blocks
template< typename T > inline udd<T> operator*(const udd<T>& a, const udd<T>& b)
{
  udd<T> c(a);
  c*=b;
  return c;
}

// divide two udd blocks
template< typename T > inline udd<T> operator/(const udd<T>& a, const udd<T>& b)
{
  udd<T> c(a);
  c/=b;
  return c;
}

// add scalar with udd block
template< typename T > inline udd<T> operator+(const T& s, const udd<T>& a)
{
  udd<T> c(a);
  c+=s;
  return c;
}

template< typename T > inline udd<T> operator+(const udd<T>&a, const T& s)
{
  udd<T> c(a);
  c+=s;
  return c;
}

// subtract scalar with udd block
template< typename T > inline udd<T> operator-(const T& s, const udd<T>& a)
{
  udd<T> c(a);
  c-=s;
  return c;
}

template< typename T > inline udd<T> operator-(const udd<T>&a, const T& s)
{
  udd<T> c(a);
  c-=s;
  return c;
}

// multiply scalar with udd block
template< typename T > inline udd<T> operator*(const T& s, const udd<T>& a)
{
  udd<T> c(a);
  c*=s;
  return c;
}

template< typename T > inline udd<T> operator*(const udd<T>& a, const T& s)
{
  udd<T> c(a);
  c*=s;
  return c;
}

// divide scalar with udd block
/*
template< typename T > inline udd<T> operator/(const T& s, const udd<T>& a)
{
  udd<T> c(a);
  c/=s;
  return c;
}
*/

template< typename T > inline udd<T> operator/(const udd<T>& a, const T& s)
{
  udd<T> c(a);
  c/=s;
  return c;
}



#endif
