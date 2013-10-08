#ifndef __CosThetaDependentInterferenceSmearer__
#define __CosThetaDependentInterferenceSmearer__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include "TFile.h"
#include "TGraphAsymmErrors.h"

class Normalization_8TeV;

// ------------------------------------------------------------------------------------
class CosThetaDependentInterferenceSmearer : public BaseGenLevelSmearer
{
public:

  CosThetaDependentInterferenceSmearer(Normalization_8TeV * norm, double &genCosTheta, std::string histFile );
  virtual ~CosThetaDependentInterferenceSmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift=0. ) const ;
  
private:
  std::string name_;
  Normalization_8TeV * norm_;
	double &genCosTheta_;
 	TFile *histFile_;
	TH1F *reweightHist_;
};

#endif
