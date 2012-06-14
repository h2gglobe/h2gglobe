#include <iostream>
#include <vector>
#include "TMath.h"

void effError(){
  
  const int nCombs=3;

  double eff[nCombs]    = {1.000,0.999,0.999};
  double effErr[nCombs] = {0.075,0.003,0.002};

  double effRes = 1.;
  double effRatSquared = 0.;

  for (int i=0; i<nCombs; i++){
    effRes*=eff[i];
    effRatSquared+=((effErr[i]/eff[i])*(effErr[i]/eff[i]));
  }

  double effErrRes = effRes*TMath::Sqrt(effRatSquared);

  cout << "Total eff: " << effRes << " +/- " << effErrRes << endl;
}
