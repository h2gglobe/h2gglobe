#ifndef MassResolution_h
#define MassResolution_h

//----------------------------------------------------------//
// Project:   MassResolution
// Author:    Matt Kenzie (matthew.william.kenzie@cern.ch)
// Modified:  25/08/2011
// Admins:    Matth Kenzie (matthew.william.kenzie@cern.ch)
//---------------------------------------------------------//

/*
  Calculates the mass resolution for a higgs to two photon event.
  Designed to work solely inside the h2gglobe framework.

  Define a MassResolution object as public data member of StatAnalysis.h
  (or equivalent) with 
    MassResolution massResFinder;
  
  In the Init method of StatAnalysis create a new object with 
    massResFinder = new MassResolution();

  Make sure to clean this in the Term method of StatAnalysis with
    delete massResFinder;

  To calculate the mass resolution event by event insert the following in 
  the Analysis method of StatAnalysis (after the smearing has been done and HiggsM andHiggs Pt have been calculated)

    massResFinder.Setup(l,diphoton_index.first,diphton_index.second,diphoton_id,higgsPt,higgsMass,eSmearPars,nR9Categories,nEtaCategories);

  and then access the mass resolution with
    
    massResFinder.getMassResolution();

  there are various other methods to get more specialised information such getAngleResolution(), getLeadResolution() etc.

  You can also use printInfo() or dumpInfo(filename) to print or write the information in Setup.

*/
  

#include <iostream>
#include <string>
#include <cassert>

#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TF1.h"
#include "math.h"

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
//#include "PhotonAnalysis.h"
#include "LoopAll.h"
#include "EnergySmearer.h"
#include "RooContainer.h"

class MassResolution {

  public:
    MassResolution();
    MassResolution(std::string);

    void Setup(LoopAll&,TLorentzVector *,TLorentzVector *,int,int,int,double,double,EnergySmearer::energySmearingParameters,int,int);

    double massResolutionCorrVtx();
    double massResolutionWrongVtx();
    double massResolutionEonly();

  private:
    double leadPhotonResolution();
    double subleadPhotonResolution();
    double getPhotonResolution(double,double,double,int,double,bool);
    
    double leadPhotonResolutionNoSmear();
    double subleadPhotonResolutionNoSmear();
    
    double angleResolution();
    double angleResolutionCorrVtx();
    double angleResolutionWrongVtx();
    
    double dzResolution();
    double dzResolutionCorrVtx();
    double dzResolutionWrongVtx();
    
   // double massResolution();
    double massResolutionAonly();
    double massResolutionEonlyNoSmear();
    double massResolutionCorrVtxNoSmear();
    double massResolutionWrongVtxNoSmear();

    //double getEffectiveSigma(TF1*, int);
    double propagateDz(double);

    double SecH(double);
    double TanH(double);
    
    void printInfo();
    void dumpInfo(std::ostream&);

    TLorentzVector *lead_p4;
    TLorentzVector *sublead_p4;
    TVector3 *lead_sc_pos;
    TVector3 *sublead_sc_pos;
    TVector3 *vertex;
    TVector3 *vtx_dxdydz;
    double lead_Eres;
    double sublead_Eres;
    double lead_r9;
    double sublead_r9;
    int lead_phoCat;
    int sublead_phoCat;
    bool lead_iDet;
    bool sublead_iDet;
    EnergySmearer::energySmearingParameters _eSmearPars;

    std::string energyCorrectionMethod;
    double dz;
    double higgsMass;

    TFile *dz_file;
    TGraph *dz_plot;

};

#endif
