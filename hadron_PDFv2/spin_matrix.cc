//  $Id$
//  $Log$
//  spin Matrix definitions
//
#include "spin_matrix.h"


SpinMatrix Cg5_matrix(){
  SpinMatrix g_one = 1.0;
  return g_one*Gamma(5) ;
}
SpinMatrix g5(){
  SpinMatrix g_one = 1.0;
  return g_one*Gamma(Ns*Ns-1) ;
}

SpinMatrix PosPar(){
  SpinMatrix g_one = 1.0;
  //spin projection matrices                                                 
  // S = (1 + gamma_4)  / 2                                                   
  //   = (1 + Gamma(8)) / 2                                                   
  return 0.5*(g_one + Gamma(8) * g_one) ;
}

SpinMatrix NegPar(){
  SpinMatrix g_one = 1.0;
  //spin projection matrices                                                  
  // S = (1 - gamma_4)  / 2                                                    
  //   = (1 - Gamma(8)) / 2                                                   
  return 0.5*(g_one - Gamma(8) * g_one) ;
}

SpinMatrix PosParSpinUp(){
  SpinMatrix g_one = 1.0;
  //spin projection matrices                                                   
  // S_proj = (1 + \Sigma_3)*(1 + gamma_4) / 2                                 
  //        = (1 + Gamma(8) - i G(3) - i G(11)) / 2                            
  return 0.5*((g_one + Gamma(8) * g_one) -
              timesI(Gamma(3) * g_one  +  Gamma(11) * g_one));
}
SpinMatrix PosParSpinDown(){
  SpinMatrix g_one = 1.0;
  //spin projection matrices                                                   
  // S_proj = (1 - \Sigma_3)*(1 + gamma_4) / 2                                 
  //        = (1 + Gamma(8) + i G(3) + i G(11)) / 2                           
  return 0.5*((g_one + Gamma(8) * g_one) +
              timesI(Gamma(3) * g_one  +  Gamma(11) * g_one));
}

SpinMatrix NegParSpinUp(){
  SpinMatrix g_one = 1.0;
  //spin projection matrices                                                  
  // S_proj = (1 + \Sigma_3)*(1 - gamma_4) / 2                                
  //        = (1 - Gamma(8) - i G(3) + i G(11)) / 2                           
  return 0.5*((g_one - Gamma(8) * g_one) -
              timesI(Gamma(3) * g_one  -  Gamma(11) * g_one));
}
SpinMatrix NegParSpinDown(){
  SpinMatrix g_one = 1.0;
  //spin projection matrices                                                  
  // S_proj = (1 - \Sigma_3)*(1 - gamma_4) / 2                                
  //        = (1 - Gamma(8) + i G(3) - i G(11)) / 2                           
  return 0.5*((g_one - Gamma(8) * g_one) +
              timesI(Gamma(3) * g_one  -  Gamma(11) * g_one));
}
