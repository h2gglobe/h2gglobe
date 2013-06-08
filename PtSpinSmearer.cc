#include "PtSpinSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

PtSpinSmearer::PtSpinSmearer(const  std::string & theFile) : efficiency_file(theFile) 
{
  name_="PtSpinSmearer";
}

PtSpinSmearer::~PtSpinSmearer()
{
}

void PtSpinSmearer::readMassPoint(int mass){
 
  ptSpinSmearersSM_[mass] = std::vector<TH1*>(3,(TH1*)0);
  ptSpinSmearersGGGRAV_[mass] = std::vector<TH1*>(3,(TH1*)0);
  ptSpinSmearersQQGRAV_[mass] = std::vector<TH1*>(3,(TH1*)0);

  // sm files
  TH1* temp = (TH1F*) thePtSpinFile_->Get( Form("smNomRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersSM_[mass][0]=(TH1F*) temp->Clone(Form("smPtMass%d",mass)); 
  ptSpinSmearersSM_[mass][0]->SetDirectory(0);

  temp = (TH1F*) thePtSpinFile_->Get( Form("smUpRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersSM_[mass][1]=(TH1F*) temp->Clone(Form("smPtMass%d_up",mass)); 
  ptSpinSmearersSM_[mass][1]->SetDirectory(0);

  temp = (TH1F*) thePtSpinFile_->Get( Form("smDownRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersSM_[mass][2]=(TH1F*) temp->Clone(Form("smPtMass%d_down",mass)); 
  ptSpinSmearersSM_[mass][2]->SetDirectory(0);

  // gg files
  temp = (TH1F*) thePtSpinFile_->Get( Form("ggNomRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersGGGRAV_[mass][0]=(TH1F*) temp->Clone(Form("ggPtMass%d",mass)); 
  ptSpinSmearersGGGRAV_[mass][0]->SetDirectory(0);

  temp = (TH1F*) thePtSpinFile_->Get( Form("ggUpRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersGGGRAV_[mass][1]=(TH1F*) temp->Clone(Form("ggPtMass%d_up",mass)); 
  ptSpinSmearersGGGRAV_[mass][1]->SetDirectory(0);

  temp = (TH1F*) thePtSpinFile_->Get( Form("ggDownRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersGGGRAV_[mass][2]=(TH1F*) temp->Clone(Form("ggPtMass%d_down",mass)); 
  ptSpinSmearersGGGRAV_[mass][2]->SetDirectory(0);

  // qq files
  temp = (TH1F*) thePtSpinFile_->Get( Form("qqNomRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersQQGRAV_[mass][0]=(TH1F*) temp->Clone(Form("qqPtMass%d",mass)); 
  ptSpinSmearersQQGRAV_[mass][0]->SetDirectory(0);

  temp = (TH1F*) thePtSpinFile_->Get( Form("qqUpRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersQQGRAV_[mass][1]=(TH1F*) temp->Clone(Form("qqPtMass%d_up",mass)); 
  ptSpinSmearersQQGRAV_[mass][1]->SetDirectory(0);

  temp = (TH1F*) thePtSpinFile_->Get( Form("qqDownRat%d",mass) );
  assert(temp!=0);    
  ptSpinSmearersQQGRAV_[mass][2]=(TH1F*) temp->Clone(Form("qqPtMass%d_down",mass)); 
  ptSpinSmearersQQGRAV_[mass][2]->SetDirectory(0);

}

bool PtSpinSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  int genMassPoint;
  string genType;

  if (sample_type==-37 || sample_type==-38 || sample_type==-39 || sample_type==-40){
    genMassPoint=125;
    genType="sm";
  }
  else if (sample_type==-77 || sample_type==-78 || sample_type==-79 || sample_type==-80){
    genMassPoint=126;
    genType="sm";
  }
  else if (sample_type==-137){
    genMassPoint=125;
    genType="gg";
  }
  else if (sample_type==-138){
    genMassPoint=125;
    genType="qq";
  }
  else if (sample_type==-177){
    genMassPoint=126;
    genType="gg";
  }
  else if (sample_type==-178){
    genMassPoint=126;
    genType="qq";
  }
  else    return true;                     // this is the case of backgrounds

  double kWeight = getWeight( p4, nPu, genMassPoint, genType, syst_shift );
  weight = (kWeight > 0) ? kWeight : 0;

  return true;
}


bool PtSpinSmearer::init() 
{
  cout << name_ << " - Opening PtSpin file"<<endl;
  thePtSpinFile_ = TFile::Open( efficiency_file.c_str() );
  assert(thePtSpinFile_!=0);

  readMassPoint(125);
  readMassPoint(126);

  thePtSpinFile_->Close();

  return true;
}

double PtSpinSmearer::getPtSpin(int genMassPoint, std::string genType, int id, double gPT  ) const {

  if(genMassPoint==125 || genMassPoint==126) {
    if (genType=="sm") {
      const TH1* tmp = ptSpinSmearersSM_.find(genMassPoint)->second[id]; 
      return tmp->GetBinContent(tmp->FindFixBin(gPT));
    }
    else if (genType=="gg"){
      const TH1* tmp = ptSpinSmearersGGGRAV_.find(genMassPoint)->second[id]; 
      return tmp->GetBinContent(tmp->FindFixBin(gPT));
    }
    else if (genType=="qq"){
      const TH1* tmp = ptSpinSmearersQQGRAV_.find(genMassPoint)->second[id]; 
      return tmp->GetBinContent(tmp->FindFixBin(gPT));
    }
    else {
      assert(0.);
      return 0.;
    }
  }
  assert(0); return 0.;
}


double PtSpinSmearer::getWeight( const TLorentzVector & p4, const int nPu, const int & genMassPoint, const std::string & genType, float syst_shift) const
{
  float gPT = p4.Pt();
  // this is consistent with samples available on Tue Jun 21 18:10:03 CEST 2011
  // bins are very fine, therefore interpolation between bins can be neglegted for now

  int    varId     = syst_shift > 0 ?   1 : 2;
  double nominal   = getPtSpin( genMassPoint, genType, 0, gPT );
  double variation = getPtSpin( genMassPoint, genType, varId, gPT );

  return ( nominal + (variation-nominal) * fabs(syst_shift) );
}
