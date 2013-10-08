#include "CosThetaDependentInterferenceSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include "TMath.h"
#include <assert.h>

#include "Macros/Normalization_8TeV.h"

CosThetaDependentInterferenceSmearer::CosThetaDependentInterferenceSmearer(Normalization_8TeV* norm, double &genCosTheta, std::string histFile) : norm_(norm), genCosTheta_(genCosTheta)
{
  name_="CosThetaDependentInterferenceSmearer";
	histFile_ = TFile::Open(histFile.c_str());
	reweightHist_ = (TH1F*)histFile_->Get("cosTheta_inteference_m125");
}

CosThetaDependentInterferenceSmearer::~CosThetaDependentInterferenceSmearer()
{
	histFile_->Close();
	delete histFile_;
}

bool CosThetaDependentInterferenceSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  int genMassPoint;

  if( sample_type >= 0 ) { return true; }
  genMassPoint = std::round(norm_->GetMass(sample_type));

  if( norm_->GetProcess(sample_type) != "ggh" && norm_->GetProcess(sample_type) != "gg_grav") {
	  return true;
  }
  if( genMassPoint > 150 ) { genMassPoint=150; } // Warning: missing k-factor
  if( genMassPoint == 100 ) { genMassPoint=105; }  // Warning: missing k-factor
 
  weight = 1.+(reweightHist_->GetBinContent(reweightHist_->FindBin(TMath::Abs(genCosTheta_)))/100.);
	//cout << "CTDepIntSmear: -- genCosTheta: " << genCosTheta_ << " wt: " << weight << endl;

	

  return true;
}


