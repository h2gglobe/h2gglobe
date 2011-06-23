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

  weight = getWeight( p4, nPu, genMassPoint, syst_shift );
  return true;
}


bool KFactorSmearer::init() 
{
  cout << name_ << " - Opening KFactor file"<<endl;
  TFile* kfacFile = TFile::Open( efficiency_file.c_str() );
  assert(kfacFile!=0);

  TH1D *temp;
  temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",110));
  assert(temp!=0);   thm110 =(TH1D*) temp->Clone(Form("Hmass%d",110)); thm110->SetDirectory(0);
  temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",120));
  assert(temp!=0);   thm120 =(TH1D*) temp->Clone(Form("Hmass%d",120)); thm120->SetDirectory(0);
  temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",130));
  assert(temp!=0);   thm130 =(TH1D*) temp->Clone(Form("Hmass%d",130)); thm130->SetDirectory(0);
  temp = (TH1D*) kfacFile->Get(Form("kfact%d_0",140));
  assert(temp!=0);   thm140 =(TH1D*) temp->Clone(Form("Hmass%d",140)); thm140->SetDirectory(0);
  
  kfacFile->Close();
  //GF get also the systematics (map?) 
  return true;
}


double KFactorSmearer::getWeight( const TLorentzVector & p4, const int nPu, const int & genMassPoint, float syst_shift) const
{
  
  float gPT = p4.Pt();

  // this is consistent with samples available on Tue Jun 21 18:10:03 CEST 2011
  // bins are very fine, therefore interpolation between bins can be neglegted for now
  if      (genMassPoint ==110 ) return thm110->GetBinContent(thm110->FindFixBin(gPT));
  else if (genMassPoint ==120 ) return thm120->GetBinContent(thm120->FindFixBin(gPT));
  else if (genMassPoint ==130 ) return thm130->GetBinContent(thm130->FindFixBin(gPT));
  else if (genMassPoint ==140 ) return thm140->GetBinContent(thm140->FindFixBin(gPT));
  else if (genMassPoint ==115 ) return (0.5*thm110->GetBinContent(thm110->FindFixBin(gPT)) +0.5*thm120->GetBinContent(thm120->FindFixBin(gPT)));
  else assert(0);

  return 1.;

}
