#ifndef __PDFWEIGHTSMEARER__
#define __PDFWEIGHTSMEARER__

#include "BaseSmearer.h"
#include <string>
#include <map>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"

class Normalization_8TeV;

// ------------------------------------------------------------------------------------
class PdfWeightSmearer : public BaseGenLevelSmearer
{
public:

  std::string efficiency_file;
  
  //   downId and upId set according to prescription of Jun 3, 2011  
  PdfWeightSmearer(const std::string &efficiency_file, Normalization_8TeV *norm, std::string downId, std::string upId) ;
  virtual ~PdfWeightSmearer();
  
  virtual const std::string & name() const { return name_; };
  
  virtual bool smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift=0. ) const ;

  void name(const std::string & x) { name_ = x; };

  void setKFName(const std::string & x) { KFName_ =x; };

  bool init();
  
 protected:
  
  double getWeight( const TLorentzVector & p4, const int nPu, float syst_shift=0.) const;
  
  std::string   name_;
  std::string   KFName_;
  Normalization_8TeV * norm_;
  TFile        *thePdfWeightFile_; 
  std::vector<TH2F*> kFactorSmearers_;
  void   readFile(std::string uId, std::string dId );
  double getPdfWeight(int genMassPoint, int id, double gPT , double gY) const;

  std::string downId, upId; 
};

#endif
