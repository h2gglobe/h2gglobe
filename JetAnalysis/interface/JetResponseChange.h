#ifndef __JetResponseChange__
#define __JetResponseChange__

#include "TH1F.h"
#include "TFile.h"
#include <string>
#include "math.h"
#include <iostream>

using namespace std;

class JetResponseChange
{
 public:
  JetResponseChange(std::string filename);
  virtual ~JetResponseChange();
  float getResponse(float eta,float lumi);
  
 private:
  TH1F* hResponse;
  
};

#endif
