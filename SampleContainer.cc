#include "SampleContainer.h"
#include <utility>
#include <iostream>

SampleContainer::~SampleContainer() 
{}

SampleContainer::SampleContainer() {
  float weight = 1;
  int itype = 0;
  int ind = 0;
  int histoplotit = 1;
  std::string filesshortnam = "";
  long long int ntot = 0;
  int nred = 0;
  float lumi = 0; 
  float xsec = 0;
  float kfactor= 1; 
  float scale = 1;
  float lumireal = 1;
}

void SampleContainer::computeWeight(float intL) {
  if(itype==0) { //this is data
    weight = 1; 
  } else {
    weight = ntot/xsec*intL;
  }
}
 
