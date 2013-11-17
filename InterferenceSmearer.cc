#include "InterferenceSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include "TMath.h"
#include <assert.h>

#include "Macros/Normalization_8TeV.h"

InterferenceSmearer::InterferenceSmearer(Normalization_8TeV* norm, double *genCosTheta, bool isConst, float correction, float error, std::string histFile) : 
	norm_(norm), genCosTheta_(genCosTheta), isConst_(isConst), correction_(correction), error_(error)
{
  name_="InterferenceSmearer";
	if (histFile != ""){
		histFile_ = TFile::Open(histFile.c_str());
		assert(histFile_);
		TH1 *temp = (TH1F*)histFile_->Get("cosTheta_interference_ggh_m125");
		assert(temp!=0);
		reweightHist_ggh_ = (TH1F*) temp->Clone("cTinterf_ggh_m125");
		reweightHist_ggh_->SetDirectory(0);
		temp = (TH1F*)histFile_->Get("cosTheta_interference_gg_grav_m125");
		assert(temp!=0);
		reweightHist_gg_grav_ = (TH1F*) temp->Clone("cTinterf_gg_grav_m125");
		reweightHist_gg_grav_->SetDirectory(0);
		histFile_->Close();
		delete histFile_;
	}
	else {
		reweightHist_ggh_ = 0;
		reweightHist_gg_grav_ = 0;
	}
}

// --- this constructor is for back compatibility ----
InterferenceSmearer::InterferenceSmearer(Normalization_8TeV* norm, float correction, float error)
{
	InterferenceSmearer(norm,NULL,true,correction,error,"");
}

InterferenceSmearer::~InterferenceSmearer()
{
	/*
	if (!isConst_) {
		histFile_->Close();
		delete histFile_;
	}
	*/
}

bool InterferenceSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  int genMassPoint;

  if( sample_type >= 0 ) { return true; }
  genMassPoint = std::round(norm_->GetMass(sample_type));

  if( norm_->GetProcess(sample_type) != "ggh" && norm_->GetProcess(sample_type) != "gg_grav") {
	  return true;
  }
  if( genMassPoint > 150 ) { genMassPoint=150; } // Warning: missing k-factor
  if( genMassPoint == 100 ) { genMassPoint=105; }  // Warning: missing k-factor

	//cout << "IntSmear: " << isConst_ << " ";
	// const interenference
	if (isConst_) {
		weight = 1.-(correction_+syst_shift*error_);
	}
	// cos theta dependent interference
	else {
		TH1 *reweightHist = 0;
		if (norm_->GetProcess(sample_type) == "ggh") reweightHist = reweightHist_ggh_;
		else if (norm_->GetProcess(sample_type) == "gg_grav") reweightHist = reweightHist_gg_grav_;
		else assert(0);
		assert(reweightHist);
		assert(genCosTheta_);
		//cout << "genCosTheta: " << *genCosTheta_ << " ";
		weight = 1.+(reweightHist->GetBinContent(reweightHist->FindBin(TMath::Abs(*genCosTheta_)))/100.);
	}
	//cout << "wt: " << weight << endl;

	

  return true;
}


