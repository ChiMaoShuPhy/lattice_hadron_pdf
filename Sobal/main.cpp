// Copyright (c) 2012 Leonhard Gruenschloss (leonhard@gruenschloss.org)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "sobol.h"
#include "sobol.cpp"
#include <cmath>
#include <time.h>
#include <fstream>
#include <iostream>

#include <vector>

#define PASS(ln) cout << "====== PASS  " << ln << " ======" << endl


using namespace std;

int nrow[4] ={4,4,4,32}; 

struct int4 {int int_arry[4];};


void sobol_src_gnrtr(int srnd_seed, int lngth, int* nrow, 
vector< int4 > &src_vct )
{
  srand(srnd_seed);
  int seed = rand()%(4*4*4*32);
  cout<<"seed = "<<seed<<endl;
  int skip = 0;
  int4 src_tmp;

   src_tmp.int_arry[0]=0; src_tmp.int_arry[1]=0;
   src_tmp.int_arry[2]=0; src_tmp.int_arry[3]=0;

  for (int i = 0; i < lngth+seed; ++i)
    {      
        // Print a few dimensions of each point.

        for (unsigned d = 0; d < 4; ++d)
        {
           const int s = (int) (nrow[d]*sobol::sample(i, d)); 
           if(skip < seed)
           { ; }
           else
           {
            src_tmp.int_arry[d] = s;                     
           }            
        }      
          if(skip < seed)
           { ; }
           else
           {
            src_vct.push_back(src_tmp);
           }
          skip++; 
       
   }
}

int main()
{
  vector<int4> src_pos;
  sobol_src_gnrtr(4, 64, nrow, src_pos);
  cout<<src_pos.size()<<endl;
  
  ofstream sobol_seqf;
    sobol_seqf.open("sobol.dat");

  for (int iter = 0; iter < src_pos.size(); iter++) 
     {  
        cout<<"No."<<iter<<":     "
            << src_pos[iter].int_arry[0]<<" "<<src_pos[iter].int_arry[1]<<" "
             << src_pos[iter].int_arry[2]<<" "<<src_pos[iter].int_arry[3]<<endl;


        sobol_seqf
              << src_pos[iter].int_arry[0]<<" "
              <<src_pos[iter].int_arry[1]<<" "
              << src_pos[iter].int_arry[2]<<" "
              <<src_pos[iter].int_arry[3]<<"\n";
     }

     sobol_seqf.close();


    
}


/*
int main(int, char**)
{

 cout<<sobol::sample(3, 2)<<endl;

 return 0;

}
*/

/*
int main(int, char**)
{

    srand(3);
    int seed = rand();
    // Iterate over points.    
    ofstream sobol_seqf;
    sobol_seqf.open("sobol.dat");

   int skip = 0;
   struct int4 {int int_arry[4];};
   vector< int4 > src_vct;
   int4 src_tmp;
   
   src_tmp.int_arry[0]=0;
   src_tmp.int_arry[1]=0;
   src_tmp.int_arry[2]=0;
   src_tmp.int_arry[3]=0;

    for (unsigned long long i = 0; i < 7+seed; ++i)
    {      
        // Print a few dimensions of each point.


        for (unsigned d = 0; d < 4; ++d)
        {
           const int s = (int) (nrow[d]*sobol::sample(i, d)); 

           if(skip < seed)
           { ; }
           else
           {
            cout << s<< " ";
          //  src_tmp[d] = s; 
            src_tmp.int_arry[d] = s;
            sobol_seqf<<s<<" ";
           }            
        }
           if(skip < seed)
           { ; }
           else
           {
            cout<<"\n";
            sobol_seqf<<"\n";
            src_vct.push_back(src_tmp);
           }

          skip++; 
    }

    sobol_seqf.close();
     cout<< "========================"<<endl;
     for (int iter = 0; iter < src_vct.size(); iter++) 
     {
       cout<<"No."<<iter<<":     "<< src_vct[iter].int_arry[0]<<" "<<src_vct[iter].int_arry[1]<<" "
             << src_vct[iter].int_arry[2]<<" "<<src_vct[iter].int_arry[3]<<endl;
     }

    return 0;
}

*/