#include "PdfWeightSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include "TMath.h"
#include <assert.h>
#include <algorithm>

#include "Macros/Normalization_8TeV.h"

PdfWeightSmearer::PdfWeightSmearer(const  std::string & theFile,  Normalization_8TeV * norm, std::string dId , std::string uId ) : 
	efficiency_file(theFile), norm_(norm), downId(dId), upId(uId)
{
  name_="PdfWeightSmearer";
}

PdfWeightSmearer::~PdfWeightSmearer()
{
}

void PdfWeightSmearer::readFile(std::string uId, std::string dId ){
  
  kFactorSmearers_.resize(3);

  TH2F* temp = (TH2F*) thePdfWeightFile_->Get("GF_cent");
  assert(temp!=0);    kFactorSmearers_[0]=(TH2F*) temp->Clone(("Hmasscent")); kFactorSmearers_[0]->SetDirectory(0);

  temp = (TH2F*) thePdfWeightFile_->Get( Form("GF_%s",uId.c_str()) );
  assert(temp!=0);    kFactorSmearers_[1]=(TH2F*) temp->Clone(("Hmass_up")); kFactorSmearers_[1]->SetDirectory(0);
  if (uId=="up"){
	// up seems to be a factor of 10 larger than it should!
	kFactorSmearers_[1]->Scale(0.1);
//	kFactorSmearers_[1]->Scale(kFactorSmearers_[0]->Integral()/kFactorSmearers_[1]->Integral());
  }

  temp = (TH2F*) thePdfWeightFile_->Get( Form("GF_%s",dId.c_str()) );
  assert(temp!=0);    kFactorSmearers_[2]=(TH2F*) temp->Clone(("Hmass_down")); kFactorSmearers_[2]->SetDirectory(0);
  if (dId=="up"){
	kFactorSmearers_[2]->Scale(0.1);
//	kFactorSmearers_[2]->Scale(kFactorSmearers_[0]->Integral()/kFactorSmearers_[2]->Integral());
  }

}

bool PdfWeightSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  if( sample_type >= 0 ) { return true; }
  int genMassPoint = std::round(norm_->GetMass(sample_type));
  
  if( norm_->GetProcess(sample_type) != "ggh" ) {
    return true;
  }
  if( genMassPoint > 150 ) { genMassPoint=150; } // Warning: missing k-factor
  if( genMassPoint == 100 ) { genMassPoint=105; }  // Warning: missing k-factor
  
  assert( genMassPoint % 5 == 0 );
  
  double kWeight = getWeight( p4, nPu, syst_shift );
  weight = (kWeight > 0) ? kWeight : 0;
  return true;
}


bool PdfWeightSmearer::init() 
{
  cout << name_ << " - Opening PdfWeight file"<<endl;
  thePdfWeightFile_ = TFile::Open( efficiency_file.c_str() );
  assert(thePdfWeightFile_!=0);

  readFile(upId, downId);
  thePdfWeightFile_->Close();

  return true;
}

double PdfWeightSmearer::getPdfWeight(int genMassPoint, int id, double gPT , double gY ) const 
{
    const TH2F* tmp = kFactorSmearers_[id]; 
    return tmp->GetBinContent(tmp->FindFixBin(gY,gPT));
    assert(0); return 0.;
}


double PdfWeightSmearer::getWeight( const TLorentzVector & p4, const int nPu, float syst_shift) const
{
  float gPT = p4.Pt();
  float gY  = fabs(p4.Rapidity());

  int    varId     = syst_shift > 0 ?   1 : 2;
  double nominal   = getPdfWeight( 0, 0, gPT, gY );
  double variation = getPdfWeight( 0, varId, gPT, gY );

  return ( max( 1. + syst_shift * (1.-variation/nominal), 0.) ) ;

}
