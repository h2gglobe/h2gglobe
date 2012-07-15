#include "PdfWeightSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include "TMath.h"
#include <assert.h>

PdfWeightSmearer::PdfWeightSmearer(const  std::string & theFile,  std::string dId , std::string uId ) : efficiency_file(theFile) , downId(dId), upId(uId)
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
  // There is a problem with the "up" histogram, its 10x as big as it should be !"
  assert(temp!=0);    kFactorSmearers_[1]=(TH2F*) temp->Clone(("Hmass_up")); kFactorSmearers_[1]->SetDirectory(0);
  if (uId=="up"){
	std::cout << "PdfWeightSmearer -- Up Histogram seems to be larger by a factor of 10 than it should be !, scaling it down!!" <<std::endl;
	kFactorSmearers_[1]->Scale(0.1);
  }

  temp = (TH2F*) thePdfWeightFile_->Get( Form("GF_%s",dId.c_str()) );
  assert(temp!=0);    kFactorSmearers_[2]=(TH2F*) temp->Clone(("Hmass_down")); kFactorSmearers_[2]->SetDirectory(0);

}

bool PdfWeightSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  // Check for GGH sample type
  if (  (sample_type == -1 )
   ||  (sample_type == -5 )
   ||  (sample_type == -9 )
   ||  (sample_type == -13)
   ||  (sample_type == -17)
   ||  (sample_type == -21)
   ||  (sample_type == -25)
   ||  (sample_type == -29)
   ||  (sample_type == -33)
   ||  (sample_type == -37)
   ||  (sample_type == -41)
   ||  (sample_type == -45)
   ||  (sample_type == -49)
   ||  (sample_type == -53)  
   ||  (sample_type == -57)  
   ||  (sample_type == -61)
   ||  (sample_type == -65)
   ||  (sample_type == -69) )
  {
 	double kWeight = getWeight( p4, nPu, syst_shift );
  	weight = (kWeight > 0) ? kWeight : 0;
	return true;
  }

  else return true; // Not a GGH type we use so  forget it
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

double PdfWeightSmearer::getPdfWeight(int genMassPoint, int id, double gPT , double gY ) const {

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

  return ( TMath::Power(variation/nominal, syst_shift) );
}
