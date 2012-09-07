#ifndef __SMEARING__
#define __SMEARING__

#include "../BaseSmearer.h"
#include "../EnergySmearer.h"
#include "../EfficiencySmearer.h"
#include "../DiPhoEfficiencySmearer.h"
#include "../KFactorSmearer.h"
#include "../InterferenceSmearer.h"


struct Smearing {
  // Add Constructor
  EfficiencySmearer::efficiencySmearingParameters effSmearPars;
  DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;
  EnergySmearer::energySmearingParameters eScaleCorrPars;
  EnergySmearer::energySmearingParameters eResolCorrPars;

  std::vector<BaseSmearer *> photonSmearers_;
  std::vector<EnergySmearer *> eScaleSmearers_;
  std::vector<EnergySmearer *> eResolSmearers_;
  std::vector<BaseDiPhotonSmearer *> diPhotonSmearers_;
  std::vector<BaseGenLevelSmearer *> genLevelSmearers_;
  std::vector<BaseSmearer *> systPhotonSmearers_;
  std::vector<BaseDiPhotonSmearer *> systDiPhotonSmearers_;
  std::vector<BaseGenLevelSmearer *> systGenLevelSmearers_;

  EnergySmearer *eScaleSmearer, *eScaleCorrSmearer;      // corrections for energy scale  MC
  EnergySmearer *eCorrSmearer;      // corrections for energy scale  MC
  EnergySmearer *eResolSmearer, *eResolCorrSmearer;
  EfficiencySmearer* idEffSmearer, *r9Smearer;
  DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer;
  
  KFactorSmearer * kFactorSmearer;
  InterferenceSmearer * interferenceSmearer;

};
#endif
