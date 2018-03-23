#include <iostream>
#include <string>
#include "chromabase.h"

using namespace Chroma ;

int main() {
  SpinMatrix gamm;
  SpinMatrix g_one = 1.0;
  int i,j,k;

  for (i=0;i<16;i++) {
    cout << "Gamma[" << i << "] :: " << endl;
    gamm = g_one * Gamma(i);
    for (j=0;j<4;j++) {
      for(k=0;k<4;k++)
	cout << peekSpin(gamm,j,k) << "\t";
      cout << endl;
    }
    cout << "\n\n" << endl;
  }

  return 0;
};
