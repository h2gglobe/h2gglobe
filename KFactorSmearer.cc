#include "KFactorSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

KFactorSmearer::KFactorSmearer(const  std::string & theFile,  unsigned int dId , unsigned int uId ) : efficiency_file(theFile) , downId(dId), upId(uId)
{
  name_="KFactorSmearer";
}

KFactorSmearer::~KFactorSmearer()
{
}

void KFactorSmearer::readMassPoint(int mass, int uId, int dId ){
  
  kFactorSmearers_[mass] = std::vector<TH1*>(3,(TH1*)0);

  TH1* temp = (TH1D*) theKFactorFile_->Get( Form("kfact%d_0",mass) );
  assert(temp!=0);    kFactorSmearers_[mass][0]=(TH1D*) temp->Clone(Form("Hmass%d",mass)); kFactorSmearers_[mass][0]->SetDirectory(0);

  temp = (TH1D*) theKFactorFile_->Get( Form("kfact%d_%d",mass,uId) );
  assert(temp!=0);    kFactorSmearers_[mass][1]=(TH1D*) temp->Clone(Form("Hmass%_up",mass)); kFactorSmearers_[mass][1]->SetDirectory(0);

  temp = (TH1D*) theKFactorFile_->Get( Form("kfact%d_%d",mass,dId) );
  assert(temp!=0);    kFactorSmearers_[mass][2]=(TH1D*) temp->Clone(Form("Hmass%_down",mass)); kFactorSmearers_[mass][2]->SetDirectory(0);

}

bool KFactorSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  int genMassPoint;

  if      (sample_type == -1 ) genMassPoint=90;
  else if (sample_type == -5 ) genMassPoint=95;
  else if (sample_type == -9 ) genMassPoint=100;
  else if (sample_type == -13) genMassPoint=105;
  else if (sample_type == -17) genMassPoint=110;
  else if (sample_type == -21) genMassPoint=115;
  else if (sample_type == -25) genMassPoint=120;
  else if (sample_type == -29) genMassPoint=130;
  else if (sample_type == -33) genMassPoint=140;
  else if (sample_type <=-37) assert(0);   // this is the case of non-existing sample
  else    return true;                     // this is the case of backgrounds

  double kWeight = getWeight( p4, nPu, genMassPoint, syst_shift );
  weight = (kWeight > 0) ? kWeight : 0;

  return true;
}


bool KFactorSmearer::init() 
{
  cout << name_ << " - Opening KFactor file"<<endl;
  theKFactorFile_ = TFile::Open( efficiency_file.c_str() );
  assert(theKFactorFile_!=0);

  readMassPoint(110, upId, downId);
  readMassPoint(120, upId, downId);
  readMassPoint(130, upId, downId);
  readMassPoint(140, upId, downId);

  theKFactorFile_->Close();

  return true;
}

double KFactorSmearer::getKFactor(int genMassPoint, int id, double gPT  ) const {

  if(genMassPoint==110 || genMassPoint==120 || genMassPoint==130  || genMassPoint==140 ) {
    const TH1* tmp = kFactorSmearers_.find(genMassPoint)->second[id]; 
    return tmp->GetBinContent(tmp->FindFixBin(gPT));
  }
  assert(0); return 0.;
}


double KFactorSmearer::getWeight( const TLorentzVector & p4, const int nPu, const int & genMassPoint, float syst_shift) const
{
  float gPT = p4.Pt();
  // this is consistent with samples available on Tue Jun 21 18:10:03 CEST 2011
  // bins are very fine, therefore interpolation between bins can be neglegted for now

  int    varId     = syst_shift > 0 ?   1 : 2;
  double nominal   = genMassPoint == 115 ?   (0.5*getKFactor(110, 0, gPT )+0.5*getKFactor(120, 0, gPT ))         : getKFactor( genMassPoint, 0, gPT );
  double variation = genMassPoint == 115 ?   (0.5*getKFactor(110, varId, gPT )+0.5*getKFactor(120, varId, gPT )) : getKFactor( genMassPoint, varId, gPT );

  return ( nominal + (variation-nominal) * fabs(syst_shift) );
}
