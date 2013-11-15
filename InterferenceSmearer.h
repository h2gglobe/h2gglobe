#ifndef __InterferenceSmearer__
#define __InterferenceSmearer__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include "TFile.h"
#include "TGraphAsymmErrors.h"

class Normalization_8TeV;

// ------------------------------------------------------------------------------------
class InterferenceSmearer : public BaseGenLevelSmearer
{
public:

  InterferenceSmearer(Normalization_8TeV * norm, double *genCosTheta, bool isConst, float correction, float error, std::string histFile );
	InterferenceSmearer(Normalization_8TeV * norm, float correction, float error);
  virtual ~InterferenceSmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift=0. ) const ;
  
private:
  std::string name_;
  Normalization_8TeV * norm_;
	double *genCosTheta_;
	bool isConst_;
	float correction_;
	float error_;
 	TFile *histFile_;
	TH1 *reweightHist_ggh_;
	TH1 *reweightHist_gg_grav_;
};

#endif
