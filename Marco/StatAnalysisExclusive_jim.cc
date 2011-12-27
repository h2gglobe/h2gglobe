#include "../interface/StatAnalysisExclusive.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
StatAnalysisExclusive::StatAnalysisExclusive()  : 
    name_("StatAnalysisExclusive"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
}

// ----------------------------------------------------------------------------------------------------
StatAnalysisExclusive::~StatAnalysisExclusive() 
{
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::Term(LoopAll& l) 
{

    std::string outputfilename = (std::string) l.histFileName;
    // Make Fits to the data-sets and systematic sets
    l.rooContainer->FitToData("data_pol_model","data_mass");  // Fit to full range of dataset
  
//    l.rooContainer->WriteSpecificCategoryDataCards(outputfilename,"data_mass","sig_mass","data_pol_model");
//    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","data_pol_model");
    // mode 0 as above, 1 if want to bin in sub range from fit,

    // Write the data-card for the Combinations Code, needs the output filename, makes binned analysis DataCard
    // Assumes the signal datasets will be called signal_name+"_mXXX"
//    l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","data_pol_model","data_mass",1,50,1); // 1 means systematics from the fit effect only the backgroundi. last digit mode = 1 means this is an internal constraint fit 
//    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","bkg_mass_rebinned");

    eventListText.close();
    eventListTextVBF.close();

    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;

    //MMMMMMM if(l.GetCutValue("optree")) 
      {
	ll->hfilereal->cd();
	optree->Write(0,TObject::kWriteDelete);
      }

//	kfacFile->Close();
//	PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::Init(LoopAll& l) 
{
  ll=&l;





  // shift and smear values
  //r9cat+nr9*etaCat
  float smear_nov14_temp[] = {0.99e-2, 1.00e-2, 1.57e-2, 2.17e-2, 3.30e-2, 3.26e-2, 3.78e-2, 3.31e-2}; 
  smear_nov14=arrayToVector(8,smear_nov14_temp);
  
  float smearErr_nov14_temp[] = {0.22e-2, 0.24e-2, 0.60e-2, 0.59e-2, 0.91e-2, 0.30e-2, 0.34e-2, 0.52e-2};  
  /*
  r9cat+nr9*etaCat
  category      minEta  MaxEta   mr9    Mr9      mrun    Mrun   smear   smearErr 
  EBlowEtaGold  0.0 1.0 0.94  999.00  -999999 999999  0.99e-2 0.22e-2
  EBlowEtaBad 0.0 1.0 -999.00 0.94  -999999 999999  1.00e-2 0.24e-2
  EBhighEtaGold 1.0 1.5 0.94  999.00  -999999 999999  1.57e-2 0.60e-2
  EBhighEtaBad  1.0 1.5 -999.00 0.94  -999999 999999  2.17e-2 0.59e-2
  EElowEtaGold  1.5 2.0 0.94  999.00  -999999 999999  3.30e-2 0.91e-2
  EElowEtaBad 1.5 2.0 -999.00 0.94  -999999 999999  3.26e-2 0.30e-2
  EEhighEtaGold 2.0 3.0 0.94  999.00  -999999 999999  3.78e-2 0.34e-2
  EEhighEtaBad  2.0 3.0 -999.00 0.94  -999999 999999  3.31e-2 0.52e-2
  */
  // runCat+nRunCat*r9cat +nRunCat*nr9Cat*etaCat

  /*
  float shift_nov14_temp[]= { 0.0028  , 0.0012  , 0.0043  , -0.0034 , -0.0038 , -0.0052, 
                    -0.0014 , -0.0030 , 0.0001  , -0.0075 , -0.0079 , -0.0094, 
                    -0.0039 , -0.0081 , -0.0039 , -0.0118 , -0.0129 , -0.0150,
                    -0.0150 , -0.0191 , -0.0150 , -0.0228 , -0.0239 , -0.0260,
                    -0.0005 , 0.0062  , 0.0048  , 0.0168  , 0.0257  , 0.0353 ,
                    -0.0025 , 0.0041  , 0.0026  , 0.0147  , 0.0236  , 0.0331 ,
                    0.0010  , 0.0113  , 0.0062  , 0.0019  , 0.0049  , 0.0009 ,
                    -0.0050 , 0.0052  , 0.0001  , -0.0041 , -0.0011 , -0.0050};
  
  shift_nov14=arrayToVector(48,shift_nov14_temp);
  */









  //if(PADEBUG) 
	cout << "InitRealStatAnalysisExclusive START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;
    
    std::string outputfilename = (std::string) l.histFileName;

    cout<<"MMMMMMMM "<<outputfilename.c_str()<<endl;

    eventListText.open(Form("%s_ascii_events.txt",outputfilename.c_str()));
    eventListTextVBF.open(Form("%s_ascii_events_vbf.txt",outputfilename.c_str()));
    FillSignalLabelMap();
    //
    // These parameters are set in the configuration file
    std::cout
	<< "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "StatAnalysis " << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "leadEtCut "<< leadEtCut << "\n"
	<< "subleadEtCut "<< subleadEtCut << "\n"
	<< "doTriggerSelection "<< doTriggerSelection << "\n"
	<< "nEtaCategories "<< nEtaCategories << "\n"
	<< "nR9Categories "<< nR9Categories << "\n"		
	<< "nPtCategories "<< nPtCategories << "\n"		
	<< "doEscaleSyst "<< doEscaleSyst << "\n"
	<< "doEresolSyst "<< doEresolSyst << "\n"
	<< "doEcorrectionSyst "<< doEcorrectionSyst << "\n"
	<< "efficiencyFile " << efficiencyFile << "\n"
	<< "doPhotonIdEffSyst "<< doPhotonIdEffSyst << "\n"
	<< "doR9Syst "<< doR9Syst << "\n"
	<< "doVtxEffSyst "<< doVtxEffSyst << "\n"
	<< "doTriggerEffSyst "<< doTriggerEffSyst << "\n"
	<< "doKFactorSyst "<< doKFactorSyst << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< std::endl;

    // avoid recalculated the CIC ID every time
    //l.runCiC = reRunCiC;
    // call the base class initializer
    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
	l.histoContainer[ind].setScale(1.);
    }
    
    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nCategories_ *= nR9Categories;
    if( nPtCategories != 0 ) nCategories_ *= nPtCategories;

    // CP

    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    
    int nVBFCategories  = ((int)includeVBF)*nVBFEtaCategories;
    int nVHadCategories = (int)includeVHad;
   
    nCategories_+=(nVBFCategories+nVHadCategories);

    l.SetCutVariables("All_phoet1",        &myAll_phoet1);
    l.SetCutVariables("All_phoet2",        &myAll_phoet2);
    l.SetCutVariables("All_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("All_phoetom2",      &myAll_phoetom2);

    l.SetCutVariables("AllLeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("AllSubJPt",        &myAllSubJPt);
    l.SetCutVariables("AllLeadJEta1",      &myAllLeadJEta);
    l.SetCutVariables("AllLeadJPt1",       &myAllLeadJPt);
    l.SetCutVariables("AllSubJPt1",        &myAllSubJPt);
    l.SetCutVariables("AllLeadJEta",      &myAllLeadJEta);
    l.SetCutVariables("AllSubJEta",       &myAllSubJEta);
    l.SetCutVariables("All_Mjj",          &myAll_Mjj);
    l.SetCutVariables("All_dEta",         &myAlldEta);
    l.SetCutVariables("All_Zep",          &myAllZep);
    l.SetCutVariables("All_dPhi",         &myAlldPhi);
    l.SetCutVariables("All_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("All_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("AllPtHiggs",       &myAllPtHiggs);

    l.SetCutVariables("All_phoet1",        &myAll_phoet1);
    l.SetCutVariables("All_phoet2",        &myAll_phoet2);
    l.SetCutVariables("All_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("All_phoetom2",      &myAll_phoetom2);

    l.SetCutVariables("All_phoet1",        &myAll_phoet1);
    l.SetCutVariables("All_phoet2",        &myAll_phoet2);
    l.SetCutVariables("All_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("All_phoetom2",      &myAll_phoetom2);

    l.SetCutVariables("All2LeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("All2SubJPt",        &myAllSubJPt);
    l.SetCutVariables("All2LeadJEta1",      &myAllLeadJEta);
    l.SetCutVariables("All2LeadJPt1",       &myAllLeadJPt);
    l.SetCutVariables("All2SubJPt1",        &myAllSubJPt);
    l.SetCutVariables("All2LeadJEta",      &myAllLeadJEta);
    l.SetCutVariables("All2SubJEta",       &myAllSubJEta);
    l.SetCutVariables("All2_Mjj",          &myAll_Mjj);
    l.SetCutVariables("All2_dEta",         &myAlldEta);
    l.SetCutVariables("All2_Zep",          &myAllZep);
    l.SetCutVariables("All2_dPhi",         &myAlldPhi);
    l.SetCutVariables("All2_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("All2_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("All2PtHiggs",       &myAllPtHiggs);

    l.SetCutVariables("All2_phoet1",        &myAll_phoet1);
    l.SetCutVariables("All2_phoet2",        &myAll_phoet2);
    l.SetCutVariables("All2_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("All2_phoetom2",      &myAll_phoetom2);


    l.SetCutVariables("Incl_Mgg0",         &myInclusive_Mgg);
    l.SetCutVariables("Incl_Mgg2",         &myInclusive_Mgg);
    l.SetCutVariables("Incl_Mgg4",         &myInclusive_Mgg);
    l.SetCutVariables("Incl_Mgg10",        &myInclusive_Mgg);
    l.SetCutVariables("InclPtHiggs",       &myInclusivePtHiggs);
    l.SetCutVariables("All_Mgg4cat",       &myInclusive_Mgg);

    l.SetCutVariables("VBF_Mgg2cat",         &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg4cat",         &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg4cat_incl",         &myAll_Mgg);

    /*
    l.SetCutVariables("VBF_phoet1",        &myAll_phoet1);
    l.SetCutVariables("VBF_phoet2",        &myAll_phoet2);
    l.SetCutVariables("VBF_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("VBF_phoetom2",      &myAll_phoetom2);
    l.SetCutVariables("VBF_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("VBFLeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("VBFSubJPt",        &myAllSubJPt);
    l.SetCutVariables("VBFLeadJEta",       &myAllLeadJEta);
    l.SetCutVariables("VBFSubJEta",        &myAllSubJEta);
    l.SetCutVariables("VBF_dEta",         &myAlldEta);
    l.SetCutVariables("VBF_Zep",          &myAllZep);
    l.SetCutVariables("VBF_Mjj",          &myAll_Mjj);
    l.SetCutVariables("VBF_dPhi",         &myAlldPhi);
    l.SetCutVariables("VBF_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg4",         &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg10",        &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg100160",    &myAll_Mgg);
    l.SetCutVariables("VBF_Mggfin",       &myAll_Mgg);
    l.SetCutVariables("VBFPtHiggs1",       &myAllPtHiggs);
    l.SetCutVariables("VBFPtHiggs2",       &myAllPtHiggs);
    l.SetCutVariables("VBFPtHiggs3",       &myAllPtHiggs);
    l.SetCutVariables("VBFPtHiggs4",       &myAllPtHiggs);
*/
    l.SetCutVariables("VBF_nvtx",        &myAll_nvtx);
    l.SetCutVariables("VBF_nvtx1",        &myAll_nvtx);
    l.SetCutVariables("VBFsr_nvtx",        &myAll_nvtx);
    l.SetCutVariables("VBFsr_nvtx1",        &myAll_nvtx);


    l.SetCutVariables("VBF_phoet1",        &myAll_phoet1);
    l.SetCutVariables("VBF_phoet2",        &myAll_phoet2);
    l.SetCutVariables("VBF_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("VBF_phoetom2",      &myAll_phoetom2);
    //l.SetCutVariables("VBF_Mgg2cat",       &myAll_Mgg);
    //l.SetCutVariables("VBF_Mgg4cat",       &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("VBFLeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("VBFSubJPt",        &myAllSubJPt);
    l.SetCutVariables("VBFLeadJEta",       &myAllLeadJEta);
    l.SetCutVariables("VBFSubJEta",        &myAllSubJEta);
    l.SetCutVariables("VBFLeadJPt1",       &myAllLeadJPt);
    l.SetCutVariables("VBFSubJPt1",        &myAllSubJPt);
    l.SetCutVariables("VBFLeadJPt2",       &myAllLeadJPt);
    l.SetCutVariables("VBFSubJPt2",        &myAllSubJPt);
    l.SetCutVariables("VBFLeadJPt3",       &myAllLeadJPt);
    l.SetCutVariables("VBFSubJPt3",        &myAllSubJPt);
    l.SetCutVariables("VBF_dEta",         &myAlldEta);
    l.SetCutVariables("VBF_Zep",          &myAllZep);
    l.SetCutVariables("VBF_Mjj",          &myAll_Mjj);
    l.SetCutVariables("VBF_dPhi",         &myAlldPhi);
    l.SetCutVariables("VBF_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg4",         &myAll_Mgg);
    l.SetCutVariables("VBF_Mgg10",        &myAll_Mgg);
    //l.SetCutVariables("VBF_Mgg100160",    &myAll_Mgg);
    //l.SetCutVariables("VBF_Mggfin",       &myAll_Mgg);
    l.SetCutVariables("VBFPtHiggs1",       &myAllPtHiggs);
    l.SetCutVariables("VBFPtHiggs2",       &myAllPtHiggs);
    l.SetCutVariables("VBFPtHiggs3",       &myAllPtHiggs);
    l.SetCutVariables("VBFPtHiggs4",       &myAllPtHiggs);

    l.SetCutVariables("VBF30_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("VBF30_phoet1",        &myAll_phoet1);
    l.SetCutVariables("VBF30_phoet2",        &myAll_phoet2);
    l.SetCutVariables("VBF30_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("VBF30_phoetom2",      &myAll_phoetom2);
    l.SetCutVariables("VBF30LeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("VBF30SubJPt",        &myAllSubJPt);
    l.SetCutVariables("VBF30LeadJEta",       &myAllLeadJEta);
    l.SetCutVariables("VBF30SubJEta",        &myAllSubJEta);
    l.SetCutVariables("VBF30_dEta",         &myAlldEta);
    l.SetCutVariables("VBF30_Zep",          &myAllZep);
    l.SetCutVariables("VBF30_Mjj",          &myAll_Mjj);
    l.SetCutVariables("VBF30_dPhi",         &myAlldPhi);
    l.SetCutVariables("VBF30_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("VBF30_Mgg10",        &myAll_Mgg);
    l.SetCutVariables("VBF30PtHiggs1",       &myAllPtHiggs);


    l.SetCutVariables("VBFsr_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("VBFsr_phoet1",        &myAll_phoet1);
    l.SetCutVariables("VBFsr_phoet2",        &myAll_phoet2);
    l.SetCutVariables("VBFsr_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("VBFsr_phoetom2",      &myAll_phoetom2);
    l.SetCutVariables("VBFsrLeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("VBFsrSubJPt",        &myAllSubJPt);
    l.SetCutVariables("VBFsrLeadJEta",       &myAllLeadJEta);
    l.SetCutVariables("VBFsrSubJEta",        &myAllSubJEta);
    l.SetCutVariables("VBFsr_dEta",         &myAlldEta);
    l.SetCutVariables("VBFsr_Zep",          &myAllZep);
    l.SetCutVariables("VBFsr_Mjj",          &myAll_Mjj);
    l.SetCutVariables("VBFsr_dPhi",         &myAlldPhi);
    l.SetCutVariables("VBFsr_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("VBFsrPtHiggs1",       &myAllPtHiggs);

    l.SetCutVariables("VBF_AB_phoet1",        &myAll_phoet1);
    l.SetCutVariables("VBF_AB_phoet2",        &myAll_phoet2);
    l.SetCutVariables("VBF_AB_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("VBF_AB_phoetom2",      &myAll_phoetom2);
    l.SetCutVariables("VBF_AB_Mgg2cat",       &myAll_Mgg);
    l.SetCutVariables("VBF_AB_Mgg4cat",       &myAll_Mgg);
    l.SetCutVariables("VBF_AB_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("VBF_ABLeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("VBF_ABSubJPt",        &myAllSubJPt);
    l.SetCutVariables("VBF_ABLeadJEta",       &myAllLeadJEta);
    l.SetCutVariables("VBF_ABSubJEta",        &myAllSubJEta);
    l.SetCutVariables("VBF_AB_dEta",         &myAlldEta);
    l.SetCutVariables("VBF_AB_Zep",          &myAllZep);
    l.SetCutVariables("VBF_AB_Mjj",          &myAll_Mjj);
    l.SetCutVariables("VBF_AB_dPhi",         &myAlldPhi);
    l.SetCutVariables("VBF_AB_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("VBF_AB_Mgg4",         &myAll_Mgg);
    l.SetCutVariables("VBF_AB_Mgg10",        &myAll_Mgg);
    l.SetCutVariables("VBF_AB_Mgg100160",    &myAll_Mgg);
    l.SetCutVariables("VBF_AB_Mggfin",       &myAll_Mgg);
    l.SetCutVariables("VBF_ABPtHiggs4",       &myAllPtHiggs);



    l.SetCutVariables("VBF1_phoet1",        &myAll_phoet1);
    l.SetCutVariables("VBF1_phoet2",        &myAll_phoet2);
    l.SetCutVariables("VBF1_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("VBF1_phoetom2",      &myAll_phoetom2);
    l.SetCutVariables("VBF1_Mgg2cat",         &myAll_Mgg);
    l.SetCutVariables("VBF1_Mgg4cat",         &myAll_Mgg);
    l.SetCutVariables("VBF1_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("VBF1LeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("VBF1SubJPt",        &myAllSubJPt);
    l.SetCutVariables("VBF1_dEta",         &myAlldEta);
    l.SetCutVariables("VBF1_Zep",          &myAllZep);
    l.SetCutVariables("VBF1_Mjj",          &myAll_Mjj);
    l.SetCutVariables("VBF1_dPhi",         &myAlldPhi);
    l.SetCutVariables("VBF1_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("VBF1_Mgg4",         &myAll_Mgg);
    l.SetCutVariables("VBF1_Mgg10",        &myAll_Mgg);
    l.SetCutVariables("VBF1_Mgg100160",    &myAll_Mgg);
    l.SetCutVariables("VBF1_Mggfin",       &myAll_Mgg);
    l.SetCutVariables("VBF1PtHiggs1",       &myAllPtHiggs);
    l.SetCutVariables("VBF1PtHiggs2",       &myAllPtHiggs);
    l.SetCutVariables("VBF1PtHiggs3",       &myAllPtHiggs);
    l.SetCutVariables("VBF1PtHiggs4",       &myAllPtHiggs);


    /*
    l.SetCutVariables("VHad_phoet1",        &myAll_phoet1);
    l.SetCutVariables("VHad_phoet2",        &myAll_phoet2);
    l.SetCutVariables("VHad_phoetom1",      &myAll_phoetom1);
    l.SetCutVariables("VHad_phoetom2",      &myAll_phoetom2);

    l.SetCutVariables("VHad_Mgg0",        &myAll_Mgg);
    l.SetCutVariables("VHadLeadJPt",      &myAllLeadJPt);
    l.SetCutVariables("VHadSubJPt",       &myAllSubJPt);
    l.SetCutVariables("VHad_dEta",        &myAlldEta);
    l.SetCutVariables("VHad_Zep",         &myAllZep);
    l.SetCutVariables("VHad_Mjj",         &myAll_Mjj);
    l.SetCutVariables("VHad_dPhi",        &myAlldPhi);
    l.SetCutVariables("VHad_Mgg2",        &myAll_Mgg);
    l.SetCutVariables("VHad_Mgg4",        &myAll_Mgg);
    l.SetCutVariables("VHad_Mgg10",       &myAll_Mgg);
    l.SetCutVariables("VHad_Mgg100160",   &myAll_Mgg);
    l.SetCutVariables("VHad_Mggfin",      &myAll_Mgg);
    */



    /*
    l.SetCutVariables("VBFR_phoet1",       &myAll_phoet1);
    l.SetCutVariables("VBFR_phoet2",       &myAll_phoet2);
    l.SetCutVariables("VBFR_phoetom1",     &myAll_phoetom1);
    l.SetCutVariables("VBFR_phoetom2",     &myAll_phoetom2);

    l.SetCutVariables("VBFR_Mgg0",         &myAll_Mgg);
    l.SetCutVariables("VBFRLeadJPt",       &myAllLeadJPt);
    l.SetCutVariables("VBFRSubJPt",        &myAllSubJPt);
    l.SetCutVariables("VBFR_dEta",         &myAlldEta);
    l.SetCutVariables("VBFR_Zep",          &myAllZep);
    l.SetCutVariables("VBFR_Mjj",          &myAll_Mjj);
    l.SetCutVariables("VBFR_dPhi",         &myAlldPhi);
    l.SetCutVariables("VBFR_Mgg2",         &myAll_Mgg);
    l.SetCutVariables("VBFR_Mgg4",         &myAll_Mgg);
    l.SetCutVariables("VBFR_Mgg10",        &myAll_Mgg);
    l.SetCutVariables("VBFR_Mgg100160",    &myAll_Mgg);
    l.SetCutVariables("VBFR_Mggfin",       &myAll_Mgg);
    */

    /*
    l.SetCutVariables("VBF_phoet1",        &myVBFphoet1);
    l.SetCutVariables("VBF_phoet2",        &myVBFphoet2);
    l.SetCutVariables("VBF_phoetom1",      &myVBFphoetom1);
    l.SetCutVariables("VBF_phoetom2",      &myVBFphoetom2);

    l.SetCutVariables("VBF_Mgg0",         &myVBF_Mgg);
    l.SetCutVariables("VBFLeadJPt",       &myVBFLeadJPt);
    l.SetCutVariables("VBFSubJPt",        &myVBFSubJPt);
    l.SetCutVariables("VBF_dEta",         &myVBFdEta);
    l.SetCutVariables("VBF_Zep",          &myVBFZep);
    l.SetCutVariables("VBF_Mjj",          &myVBF_Mjj);
    l.SetCutVariables("VBF_dPhi",         &myVBFdPhi);
    l.SetCutVariables("VBF_Mgg2",         &myVBF_Mgg);
    l.SetCutVariables("VBF_Mgg4",         &myVBF_Mgg);
    l.SetCutVariables("VBF_Mgg10",        &myVBF_Mgg);
    l.SetCutVariables("VBF_Mgg100160",    &myVBF_Mgg);
    l.SetCutVariables("VBF_Mggfin",       &myVBF_Mgg);

    l.SetCutVariables("VHad_phoet1",        &myVHadphoet1);
    l.SetCutVariables("VHad_phoet2",        &myVHadphoet2);
    l.SetCutVariables("VHad_phoetom1",      &myVHadphoetom1);
    l.SetCutVariables("VHad_phoetom2",      &myVHadphoetom2);

    l.SetCutVariables("VHad_Mgg0",        &myVHad_Mgg);
    l.SetCutVariables("VHadLeadJPt",      &myVHadLeadJPt);
    l.SetCutVariables("VHadSubJPt",       &myVHadSubJPt);
    l.SetCutVariables("VHad_dEta",        &myVHaddEta);
    l.SetCutVariables("VHad_Zep",         &myVHadZep);
    l.SetCutVariables("VHad_Mjj",         &myVHad_Mjj);
    l.SetCutVariables("VHad_dPhi",        &myVHaddPhi);
    l.SetCutVariables("VHad_Mgg2",        &myVHad_Mgg);
    l.SetCutVariables("VHad_Mgg4",        &myVHad_Mgg);
    l.SetCutVariables("VHad_Mgg10",       &myVHad_Mgg);
    l.SetCutVariables("VHad_Mgg100160",   &myVHad_Mgg);
    l.SetCutVariables("VHad_Mggfin",      &myVHad_Mgg);

    l.SetCutVariables("VBFR_phoet1",       &myVBFRphoet1);
    l.SetCutVariables("VBFR_phoet2",       &myVBFRphoet2);
    l.SetCutVariables("VBFR_phoetom1",     &myVBFRphoetom1);
    l.SetCutVariables("VBFR_phoetom2",     &myVBFRphoetom2);

    l.SetCutVariables("VBFR_Mgg0",         &myVBF_Mgg);
    l.SetCutVariables("VBFRLeadJPt",       &myVBFLeadJPt);
    l.SetCutVariables("VBFRSubJPt",        &myVBFSubJPt);
    l.SetCutVariables("VBFR_dEta",         &myVBFdEta);
    l.SetCutVariables("VBFR_Zep",          &myVBFZep);
    l.SetCutVariables("VBFR_Mjj",          &myVBF_Mjj);
    l.SetCutVariables("VBFR_dPhi",         &myVBFdPhi);
    l.SetCutVariables("VBFR_Mgg2",         &myVBF_Mgg);
    l.SetCutVariables("VBFR_Mgg4",         &myVBF_Mgg);
    l.SetCutVariables("VBFR_Mgg10",        &myVBF_Mgg);
    l.SetCutVariables("VBFR_Mgg100160",    &myVBF_Mgg);
    l.SetCutVariables("VBFR_Mggfin",       &myVBF_Mgg);
*/

    /*
    l.SetCutVariables("VHadMLeadJPt",      &myVHadLeadJPt);
    l.SetCutVariables("VHadMSubJPt",       &myVHadSubJPt);
    l.SetCutVariables("VHadM_dEta",        &myVHaddEta);
    l.SetCutVariables("VHadM_Zep",         &myVHadZep);
    l.SetCutVariables("VHadM_Mjj",         &myVHad_Mjj);
    l.SetCutVariables("VHadM_dPhi",        &myVHaddPhi);
    l.SetCutVariables("VHadM_Mgg0",        &myVHad_Mgg);
    l.SetCutVariables("VHadM_Mgg2",        &myVHad_Mgg);
    l.SetCutVariables("VHadM_Mgg4",        &myVHad_Mgg);
    l.SetCutVariables("VHadM_Mgg10",       &myVHad_Mgg);
    l.SetCutVariables("VHadM_Mgg100160",   &myVHad_Mgg);
    l.SetCutVariables("VHadM_Mggfin",      &myVHad_Mgg);
    */

    // CP


    effSmearPars.categoryType = "2CatR9_EBEE";
    effSmearPars.n_categories = 4;
    effSmearPars.efficiency_file = efficiencyFile;

    diPhoEffSmearPars.n_categories = 8;
    diPhoEffSmearPars.efficiency_file = efficiencyFile;

    if( doEcorrectionSmear ) {
        // instance of this smearer done in PhotonAnalysis
        photonSmearers_.push_back(eCorrSmearer);
    }
    if( doEscaleSmear ) {
        photonSmearers_.push_back(eScaleSmearer);
    }
    if( doEresolSmear ) {
	// energy resolution smearing
	std::cerr << __LINE__ << std::endl; 
	eResolSmearer = new EnergySmearer( eSmearPars );
	eResolSmearer->name("E_res");
	eResolSmearer->doEnergy(false);
	eResolSmearer->scaleOrSmear(false);
	photonSmearers_.push_back(eResolSmearer);
    }
    if( doPhotonIdEffSmear ) {
	// photon ID efficiency 
	std::cerr << __LINE__ << std::endl; 
	idEffSmearer = new EfficiencySmearer( effSmearPars );
	idEffSmearer->name("idEff");
	idEffSmearer->setEffName("ratioTP");
	idEffSmearer->init();
	idEffSmearer->doPhoId(true);
	photonSmearers_.push_back(idEffSmearer);
    }
    if( doR9Smear ) {
	// R9 re-weighting
	r9Smearer = new EfficiencySmearer( effSmearPars );
	r9Smearer->name("r9Eff");
	r9Smearer->setEffName("ratioR9");
	r9Smearer->init();
	r9Smearer->doR9(true);
	photonSmearers_.push_back(r9Smearer);
    }
    if( doVtxEffSmear ) {
	// Vertex ID
	std::cerr << __LINE__ << std::endl; 
	vtxEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );   // triplicate TF1's here
	vtxEffSmearer->name("vtxEff");
	vtxEffSmearer->setEffName("ratioVertex");
	vtxEffSmearer->doVtxEff(true);
	vtxEffSmearer->init();
	diPhotonSmearers_.push_back(vtxEffSmearer);
    }
    if( doTriggerEffSmear ) {
	// trigger efficiency
	std::cerr << __LINE__ << std::endl; 
	triggerEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
	triggerEffSmearer->name("triggerEff");
	triggerEffSmearer->setEffName("effL1HLT");
	triggerEffSmearer->doVtxEff(false);
	triggerEffSmearer->init();
	diPhotonSmearers_.push_back(triggerEffSmearer);
    }
    if(doKFactorSmear) {
	// kFactor efficiency
	std::cerr << __LINE__ << std::endl; 
	kFactorSmearer = new KFactorSmearer( kfacHist );
	kFactorSmearer->name("kFactor");
	kFactorSmearer->init();
	genLevelSmearers_.push_back(kFactorSmearer);
    }

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed

    // FIXME move these params to config file
    l.rooContainer->SetNCategories(nCategories_);
    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;
    // RooContainer does not support steps different from 1 sigma
    //assert( ((float)nSystSteps) == systRange );
    if( doEcorrectionSmear && doEcorrectionSyst ) {
        // instance of this smearer done in PhotonAnalysis
        systPhotonSmearers_.push_back(eCorrSmearer);
	std::vector<std::string> sys(1,eCorrSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEscaleSmear && doEscaleSyst ) {
	systPhotonSmearers_.push_back( eScaleSmearer );
	std::vector<std::string> sys(1,eScaleSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
	systPhotonSmearers_.push_back( eResolSmearer );
	std::vector<std::string> sys(1,eResolSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonIdEffSmear && doPhotonIdEffSyst ) {
	systPhotonSmearers_.push_back( idEffSmearer );
	std::vector<std::string> sys(1,idEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doR9Smear && doR9Syst ) {
	systPhotonSmearers_.push_back( r9Smearer );
	std::vector<std::string> sys(1,r9Smearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doVtxEffSmear && doVtxEffSyst ) {
	systDiPhotonSmearers_.push_back( vtxEffSmearer );
	std::vector<std::string> sys(1,vtxEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doTriggerEffSmear && doTriggerEffSyst ) {
	systDiPhotonSmearers_.push_back( triggerEffSmearer );
	std::vector<std::string> sys(1,triggerEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if(doKFactorSmear && doKFactorSyst) {
	systGenLevelSmearers_.push_back(kFactorSmearer);
	std::vector<std::string> sys(1,kFactorSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
	
    // ----------------------------------------------------
    // ----------------------------------------------------
    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.045,1.00);
    // ----------------------------------------------------

    // Create observables for shape-analysis with ranges
    // l.rooContainer->AddObservable("mass" ,100.,150.);
    l.rooContainer->AddObservable("CMS_hgg_mass" ,massMin,massMax);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);

    // SM Model
    l.rooContainer->AddConstant("XSBR_ggh_150",0.01428);
    l.rooContainer->AddConstant("XSBR_vbf_150",0.001308);
    l.rooContainer->AddConstant("XSBR_wzh_150",0.000641);
    l.rooContainer->AddConstant("XSBR_tth_150",0.000066);
    l.rooContainer->AddConstant("XSBR_ggh_145",0.018820);
    l.rooContainer->AddConstant("XSBR_vbf_145",0.001676);
    l.rooContainer->AddConstant("XSBR_wzh_145",0.000891);
    l.rooContainer->AddConstant("XSBR_tth_145",0.000090);
    l.rooContainer->AddConstant("XSBR_ggh_140",0.0234109);
    l.rooContainer->AddConstant("XSBR_vbf_140",0.00203036);
    l.rooContainer->AddConstant("XSBR_wzh_140",0.001163597);
    l.rooContainer->AddConstant("XSBR_tth_140",0.000117189);
    l.rooContainer->AddConstant("XSBR_ggh_135",0.0278604);
    l.rooContainer->AddConstant("XSBR_vbf_135",0.002343);
    l.rooContainer->AddConstant("XSBR_wzh_135",0.001457559);
    l.rooContainer->AddConstant("XSBR_tth_135",0.000145053);
    l.rooContainer->AddConstant("XSBR_ggh_130",0.0319112);
    l.rooContainer->AddConstant("XSBR_vbf_130",0.00260804);
    l.rooContainer->AddConstant("XSBR_wzh_130",0.001759636);
    l.rooContainer->AddConstant("XSBR_tth_130",0.000173070);
    l.rooContainer->AddConstant("XSBR_ggh_125",0.0350599);
    l.rooContainer->AddConstant("XSBR_vbf_125",0.00277319);
    l.rooContainer->AddConstant("XSBR_wzh_125",0.002035123);
    l.rooContainer->AddConstant("XSBR_tth_125",0.000197718);
    l.rooContainer->AddConstant("XSBR_ggh_120",0.0374175);
    l.rooContainer->AddConstant("XSBR_vbf_120",0.00285525);
    l.rooContainer->AddConstant("XSBR_wzh_120",0.002285775);
    l.rooContainer->AddConstant("XSBR_tth_120",0.00021951);
    l.rooContainer->AddConstant("XSBR_ggh_123",0.0360696);
    l.rooContainer->AddConstant("XSBR_vbf_123",0.00281352);
    l.rooContainer->AddConstant("XSBR_wzh_123",0.00213681);
    l.rooContainer->AddConstant("XSBR_tth_123",0.00020663);
    l.rooContainer->AddConstant("XSBR_ggh_121",0.0369736);
    l.rooContainer->AddConstant("XSBR_vbf_121",0.00284082);
    l.rooContainer->AddConstant("XSBR_wzh_121",0.00223491);
    l.rooContainer->AddConstant("XSBR_tth_121",0.00021510);
    l.rooContainer->AddConstant("XSBR_ggh_115",0.0386169);
    l.rooContainer->AddConstant("XSBR_vbf_115",0.00283716);
    l.rooContainer->AddConstant("XSBR_wzh_115",0.002482089);
    l.rooContainer->AddConstant("XSBR_tth_115",0.000235578);
    l.rooContainer->AddConstant("XSBR_ggh_110",0.0390848);
    l.rooContainer->AddConstant("XSBR_vbf_110",0.00275406);
    l.rooContainer->AddConstant("XSBR_wzh_110",0.002654575);
    l.rooContainer->AddConstant("XSBR_tth_110",0.000247629);
    l.rooContainer->AddConstant("XSBR_ggh_105",0.0387684);
    l.rooContainer->AddConstant("XSBR_vbf_105",0.00262016);
    l.rooContainer->AddConstant("XSBR_wzh_105",0.002781962);
    l.rooContainer->AddConstant("XSBR_tth_105",0.000255074);

    // FF model	
    l.rooContainer->AddConstant("ff_XSBR_vbf_150",0.00259659);
    l.rooContainer->AddConstant("ff_XSBR_wzh_150",0.00127278);
    l.rooContainer->AddConstant("ff_XSBR_vbf_145",0.00387544);
    l.rooContainer->AddConstant("ff_XSBR_wzh_145",0.00205969);
    l.rooContainer->AddConstant("ff_XSBR_vbf_140",0.00565976);
    l.rooContainer->AddConstant("ff_XSBR_wzh_140",0.003243602);
    l.rooContainer->AddConstant("ff_XSBR_vbf_135",0.00825);
    l.rooContainer->AddConstant("ff_XSBR_wzh_135",0.00513225);
    l.rooContainer->AddConstant("ff_XSBR_vbf_130",0.0122324);
    l.rooContainer->AddConstant("ff_XSBR_wzh_130",0.00825316);
    l.rooContainer->AddConstant("ff_XSBR_vbf_125",0.0186494);
    l.rooContainer->AddConstant("ff_XSBR_wzh_125",0.01368598);
    l.rooContainer->AddConstant("ff_XSBR_vbf_123",0.022212);
    l.rooContainer->AddConstant("ff_XSBR_wzh_123",0.0168696);
    l.rooContainer->AddConstant("ff_XSBR_vbf_121",0.0266484);
    l.rooContainer->AddConstant("ff_XSBR_wzh_121",0.0209646);
    l.rooContainer->AddConstant("ff_XSBR_vbf_120",0.0293139);
    l.rooContainer->AddConstant("ff_XSBR_wzh_120",0.02346729);
    l.rooContainer->AddConstant("ff_XSBR_vbf_115",0.0482184);
    l.rooContainer->AddConstant("ff_XSBR_wzh_115",0.04218386);
    l.rooContainer->AddConstant("ff_XSBR_vbf_110",0.083181);
    l.rooContainer->AddConstant("ff_XSBR_wzh_110",0.08017625);
    l.rooContainer->AddConstant("ff_XSBR_vbf_105",0.151616);
    l.rooContainer->AddConstant("ff_XSBR_wzh_105",0.1609787);

    // Background modeling 
    l.rooContainer->AddRealVar("pol0",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol1",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol2",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modpol0","@0*@0","pol0");
    l.rooContainer->AddFormulaVar("modpol1","@0*@0","pol1");
    l.rooContainer->AddFormulaVar("modpol2","@0*@0","pol2");

    l.rooContainer->AddRealVar("expol0",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("expol1",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modexpol0","@0*@0","expol0");
    l.rooContainer->AddFormulaVar("modexpol1","@0*@0","expol1");
    
    // Generic PDF ok in the std analysis but excluisve channels need different models CP
    //l.rooContainer->AddGenericPdf("data_pol_model",
	  //"0","CMS_hgg_mass",data_pol_pars,73);	// >= 71 means RooBernstein of order >= 1
        

    int cats_with_std[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int cats_with_excl[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    for(int i=0; i<nEtaCategories*nR9Categories*nPtCategories; i++){
      cats_with_std[i]=1;
      cats_with_excl[i]=0;
    }


    std::vector<std::string> data_pol_pars(3,"p");	 
    data_pol_pars[0] = "modpol0";
    data_pol_pars[1] = "modpol1";
    data_pol_pars[2] = "modpol2";
    l.rooContainer->AddSpecificCategoryPdf(cats_with_std,"data_pol_model",
	  "0","CMS_hgg_mass",data_pol_pars,73);	// >= 71 means RooBernstein of order >= 1
    

    std::vector<std::string> data_excl_pars(2,"p");	 
    data_excl_pars[0] = "modexpol0";
    data_excl_pars[1] = "modexpol1";
    l.rooContainer->AddSpecificCategoryPdf(cats_with_excl, "data_pol_model",
	  "0","CMS_hgg_mass",data_excl_pars,72);	// >= 71 means RooBernstein of order >= 1

    // CP

    // -----------------------------------------------------
    // Make some data sets from the observables to fill in the event loop		  
    // Binning is for histograms (will also produce unbinned data sets)

    l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass"    ,nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass"     ,nDataBins);    	  	

    // Create Signal DataSets:
    for (int sig=105;sig<=150;sig+=5){
	// Needed to use S4 for the GGH 145 Signal which has the BUG so no 145 sample
	if (sig==145) continue;
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_ggh_mass_m%d",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_vbf_mass_m%d",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_wzh_mass_m%d",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_tth_mass_m%d",sig),nDataBins);   
 
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_ggh_mass_m%d_rv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_vbf_mass_m%d_rv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_wzh_mass_m%d_rv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_tth_mass_m%d_rv",sig),nDataBins);    
                                                                              
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_ggh_mass_m%d_wv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_vbf_mass_m%d_wv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_wzh_mass_m%d_wv",sig),nDataBins);    
    	l.rooContainer->CreateDataSet("CMS_hgg_mass",Form("sig_tth_mass_m%d_wv",sig),nDataBins);    
    }

    // Also create the 121 and 123 test points
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m121",nDataBins);   
 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m121_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m121_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m121_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m121_rv",nDataBins);    
                                                                      
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m121_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m121_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m121_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m121_wv",nDataBins);    
	
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m123",nDataBins);   
 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m123_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m123_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m123_rv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m123_rv",nDataBins);    
                                                                      
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_ggh_mass_m123_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_vbf_mass_m123_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_wzh_mass_m123_wv",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_tth_mass_m123_wv",nDataBins);    


    // Make more datasets representing Systematic Shifts of various quantities

    for (int sig=105;sig<=150;sig+=5){
	if (sig==145) continue;
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_ggh_mass_m%d",sig),-1);	
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_vbf_mass_m%d",sig),-1);	
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_wzh_mass_m%d",sig),-1);	
    	l.rooContainer->MakeSystematics("CMS_hgg_mass",Form("sig_tth_mass_m%d",sig),-1);	
    }

    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_ggh_mass_m121",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_vbf_mass_m121",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_wzh_mass_m121",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_tth_mass_m121",-1);
	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_ggh_mass_m123",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_vbf_mass_m123",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_wzh_mass_m123",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_tth_mass_m123",-1);	

    // Make sure the Map is filled
    FillSignalLabelMap();






    if(l.GetCutValue("optree")) {


    //JIM START OF THINGS
 

      photonCutsSet4=false;
      photonCutsSet6=false;
      photonCutsSet6pf=false;
      
      //DODELTARCUT=false;
      
      //GENMATCH_DELTAR_THRESHOLD = 0.1;
      
      //nduplicate=0;
      
      
      if(l.typerun==0) {
	setDiphoCuts();
	if(l.GetCutValue("bdt")) {
	  SetBDT();
	}                            // jgb 
	for(int idc=0;idc!=13;++idc) {
	  std::cout << "\ncutdiphoptom         ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutdiphoptom[idc][ic];}
	  std::cout << "\ncutfsubleadcutindex  ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutfsubleadcutindex[idc][ic];}
	  std::cout << "\ncutfleadcutindex     ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutfleadcutindex[idc][ic];}
	  std::cout << "\ncutetamax            ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutetamax[idc][ic];}
	  std::cout << "\ncutetamin            ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutetamin[idc][ic];}
	  std::cout << "\ncutsubleadptomass    ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutsubleadptomass[idc][ic];}
	  std::cout << "\ncutleadptomass       ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutleadptomass[idc][ic];}
	  std::cout << "\ncutsumptom              ";
	  for(int ic=0;ic!=4;++ic){std::cout << "\t["<<idc<<"]["<<ic<<"]: " << cutsumptom[idc][ic];}
	  std::cout << std::                    endl;
	}
	
	for(int icutlevel=0;icutlevel!=phoNCUTLEVELS;++icutlevel) {
	  std::cout << "CUT LEVEL: " << icutlevel << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_isosumoet[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_isosumoetbad[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_trkisooet[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_sieie[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_hovere[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_r9[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_drtotk_25_99[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	  for(int icat=0;icat!=4;++icat) std::cout << Ci4C_cut_lead_pixel[icutlevel][icat] << "\t";
	  std::cout << std::endl;
	}
	
	
	for(int iiii=0;iiii!=21;++iiii) {
	  allevtvtxcount[iiii]=0;
	  goodevtvtxcount1[iiii]=0;
	  goodevtvtxcount14[iiii]=0;
	  goodevtvtxcount17[iiii]=0;
	  goodevtvtxcount2[iiii]=0;
	}
	
      //END JIM
      }

    }
      
    //if(l.GetCutValue("optree")) 
    {
      HggBookOptree();
    }
    
    if(PADEBUG) 
      cout << "InitRealStatAnalysisExclusive END"<<endl;
	
    // FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::Analysis(LoopAll& l, Int_t jentry) 
{
    if(PADEBUG) 
	cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
   
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;
    float newweight = l.sampleContainer[l.current_sample_index].weight;
    double pileupWeight=1.; 



    l.FillCounter( "Processed", 1. );
    assert( weight > 0. );  
    l.FillCounter( "XSWeighted", weight );
    nevents+=1.;

    //PU reweighting
    unsigned int n_pu = l.pu_n;
    if ( cur_type !=0 && puHist != "") {
	bool hasSpecificWeight = weights.find( cur_type ) != weights.end() ; 
	if( cur_type < 0 && !hasSpecificWeight && jentry == 1 ) {
	    std::cerr  << "WARNING no pu weights specific for sample " << cur_type << std::endl;
	}
	std::vector<double> & puweights = hasSpecificWeight ? weights[ cur_type ] : weights[0]; 
	if(n_pu<puweights.size()){
	  weight *= puweights[n_pu]; //MARCO
	  pileupWeight *= puweights[n_pu];

	    sumwei+=puweights[n_pu]; 
	}    
	else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	    cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< l.itype[l.current]<<"], event will not be reweighted for pileup"<<endl;
	}
    }
    
    assert( weight >= 0. );  //marco
    l.FillCounter( "PUWeighted", weight ); //marco
    
    if( jentry % 10000 ==  0 ) {
	    std::cout << " nevents " <<  nevents << " sumpuweights " << sumwei << " ratio " << sumwei / nevents 
		      << " equiv events " << sumev << " accepted " << sumaccept << " smeared " << sumsmear << " "  
		      <<  sumaccept / sumev << " " << sumsmear / sumaccept
		      << std::endl;
    }
    // ------------------------------------------------------------
    //PT-H K-factors
    double gPT = 0;
    TLorentzVector gP4(0,0,0,0);
    if (cur_type<0){            // if background sample, gP4 remains 4vect(0)
	for (int gi=0;gi<l.gp_n;gi++){
	    if (l.gp_pdgid[gi]==25){
		gP4 = *((TLorentzVector*)l.gp_p4->At(gi));
		gPT = gP4.Pt();
		break;
	    }
	}
    }

    // ------------------------------------------------------------

    // smear all of the photons!
    std::pair<int,int> diphoton_index;
   
    // do gen-level dependent first (e.g. k-factor); only for signal
    double genLevWeight=1; 
    if(cur_type!=0){
	for(std::vector<BaseGenLevelSmearer*>::iterator si=genLevelSmearers_.begin(); si!=genLevelSmearers_.end(); si++){
	    float genWeight=1;
	    (*si)->smearEvent( genWeight,gP4, l.pu_n, cur_type, 0. );
	    if( genWeight < 0. ) {
		std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		assert(0);
	    }
	    genLevWeight*=genWeight;
	}
    }

    // Nominal smearing
    std::vector<float> smeared_pho_energy(l.pho_n,0.); 
    std::vector<float> smeared_pho_r9(l.pho_n,0.); 
    std::vector<float> smeared_pho_weight(l.pho_n,1.);

    // TEMPORARY FIX -------------------------------------------------------------------------------------------------------//
    // Scale all the r9 of the photons in the MC
    // For now we just let it use the index but we specifically Change the r9 in the branch AFTER Energy regression smearing
    // Ideally we want to pass a smeared r9 too and apply after energy corrections, currently the smeared_pho_r9 isnt used!
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    if (cur_type !=0){
      for (int ipho=0;ipho<l.pho_n;ipho++){
        double R9_rescale = (l.pho_isEB[ipho]) ? 1.0048 : 1.00492 ;
        l.pho_r9[ipho]*=R9_rescale; //commented MARCO for now, should ask
      }
    }
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
   
    for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
	std::vector<std::vector<bool> > p;
	PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
				    // *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
				    ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
				    energyCorrected[ipho],
				    l.pho_isEB[ipho], l.pho_r9[ipho],
				    l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_),
				    (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
				    );
	float pweight = 1.;
	// smear MC. But apply energy shift to data 
	if( cur_type != 0 && doMCSmearing ) { // if it's MC
	    for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
		float sweight = 1.;
		(*si)->smearPhoton(phoInfo,sweight,l.run,0.);	   
		if( sweight < 0. ) {
			std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
			assert(0);
		}
		pweight *= sweight;
	    }
	} else if( cur_type == 0 ) {          // if it's data
	    float sweight = 1.;
	    if( doEcorrectionSmear )  { 
	      eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.); 
	    }
	    eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
	    pweight *= sweight;
	}
	smeared_pho_energy[ipho] = phoInfo.energy();
	smeared_pho_r9[ipho] = phoInfo.r9();
	smeared_pho_weight[ipho] = pweight;
    }
   
    sumev += weight;

    int ccat1=-1;
    int ccat2=-1;
    int ccat3=-1;
    int ccat4=-1;



    if(l.GetCutValue("optree")) {
      //here jim's stuff for now:

      //HERE WE GET THE diphoton id and vtxind from the standard one??? //MARCO FIX
      
      int diphoton_id_jim = l.DiphotonCiCSelection(l.phoLOOSE, l.phoLOOSE, 40., 30., 4,applyPtoverM, &smeared_pho_energy[0] ); 
      if(diphoton_id_jim>=0)
	diphoton_id_jim = l.DiphotonCiCSelection(l.phoTIGHT, l.phoTIGHT, 40., 30., 4,applyPtoverM, &smeared_pho_energy[0] ); 
      
      if(diphoton_id_jim>=0) {
	TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id_jim], l.dipho_vtxind[diphoton_id_jim], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id_jim], l.dipho_vtxind[diphoton_id_jim], &smeared_pho_energy[0]);
	TLorentzVector diphoton = lead_p4+sublead_p4;
	float evweight = newweight * smeared_pho_weight[l.dipho_leadind[diphoton_id_jim]] * smeared_pho_weight[l.dipho_subleadind[diphoton_id_jim]] * genLevWeight * pileupWeight;
	
	if(fabs((float) newweight*pileupWeight-weight)/((float) newweight*pileupWeight+weight)>0.0001) cout<<"################ "<<newweight*pileupWeight<<" "<<weight<<" "<<newweight<<" "<<pileupWeight<<endl;
	
	float myweight=1.;
	if(evweight*newweight!=0) myweight=evweight/newweight;
	
	
	//check itype
	
	//NEED TO PASS THE WEIGHTS
	
	int itypepass=cur_type;
	
	SetOutputNtupleVariables(jentry, itypepass, l.dipho_leadind[diphoton_id_jim], l.dipho_subleadind[diphoton_id_jim], l.dipho_vtxind[diphoton_id_jim], diphoton.M(), &lead_p4, &sublead_p4, evweight, pileupWeight);

	t_pvtx = vtxAna_.vertexProbability((*l.vtx_std_evt_mva)[l.dipho_vtxind[diphoton_id_jim]]);
        Float_t t_dmodz = l.getDmOverDz(l.dipho_leadind[diphoton_id_jim], l.dipho_subleadind[diphoton_id_jim], &smeared_pho_energy[0]);
        float z_gg = ((TVector3*)(l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id_jim])))->Z();
        t_sigma_mz = fabs(t_dmodz)*(sqrt(pow(double(l.bs_sigmaZ), 2) + pow(double(z_gg), 2))) / diphoton.M();
        t_bsZ = l.bs_sigmaZ;

	t_leadphoidmitmva = l.photonIDMVA(l.dipho_leadind[diphoton_id_jim], l.dipho_vtxind[diphoton_id_jim], "MIT");
	t_subleadphoidmitmva = l.photonIDMVA(l.dipho_subleadind[diphoton_id_jim], l.dipho_vtxind[diphoton_id_jim], "MIT");
        t_diphomitmva = l.diphotonMVA(l.dipho_leadind[diphoton_id_jim], l.dipho_subleadind[diphoton_id_jim], l.dipho_vtxind[diphoton_id_jim], t_pvtx, diphoton.Pt(), diphoton.M(), "MIT");

	optree->Fill();
	int noptree=optree->GetEntries();
	if(noptree<100) cout<<" filling optree n="<<noptree<<endl;
	if(noptree%1000==1) cout<<" optree n="<<noptree<<endl;
	//if(noptree%10000==2) optree->Print();
	
      }
    }


    // FIXME pass smeared R9
    int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 

    int passincl=0;
    int catincl=0;
    float evweightincl=0.;
    float myweightincl=0.;
    float evweightvbf=0.;
    float myweightvbf=0.;

    if(diphoton_id>-1)
    {
      TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
      TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
      TLorentzVector diphoton = lead_p4+sublead_p4;
      myAll_Mgg =diphoton.M();
      myInclusive_Mgg = diphoton.M();
      myInclusivePtHiggs =diphoton.Pt();
      
      //should be done for both earlier
      diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
      float evweight = newweight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight * pileupWeight;

      if(fabs((float) newweight*pileupWeight-weight)/((float) newweight*pileupWeight+weight)>0.0001) cout<<"################ "<<newweight*pileupWeight<<" "<<weight<<" "<<newweight<<" "<<pileupWeight<<endl;
      //if(newweight*pileupWeight != weight) cout<<"################ "<<newweight*pileupWeight<<" "<<weight<<" "<<newweight<<" "<<pileupWeight<<endl;

      float myweight=1.;
      if(evweight*newweight!=0) myweight=evweight/newweight;
      
	if(l.ApplyCut("massrange",myAll_Mgg,0)) {
	  l.ApplyCutsFill(0,4,evweight, myweight);
	  passincl=1;
	}
      if(myInclusive_Mgg>100.&&myInclusive_Mgg<180.) {
	l.FillHist("run",0, l.run, 1.);
      }
      evweightincl=evweight;
      myweightincl=myweight;


      int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,0.,2,2,1);
      catincl=category;

      ccat1=category;

      //cout<<"MMM "<<category <<" "<<evweight<<" "<<myweight<<endl;

      if(myInclusive_Mgg>100&&myInclusive_Mgg<180) {
	//cut_All_Mgg4cat
	l.ApplyCutsFill(category,14,evweight, myweight);
      }

      //cout<<"MMM1 "<<category <<" "<<evweight<<" "<<myweight<<endl;
    }
    //cout<<"MMM1 "<<endl;

    // CP
    int diphotonVBF_id = -1;
    int diphotonVHad_id = -1;
    bool VBFevent = false;
    bool VHadevent = false;
    std::pair<int,int> highestPtJets(-1,-1);
    //if((includeVBF || includeVHad)&&l.jet_algoPF1_n>1) {
    {
      if(l.jet_algoPF1_n>1) {
	RescaleJetEnergy(l);
      }

      int applyPtoverM=0;

      leadEtVBFCut=40.;
      leadEtVHadCut=40.;
      subleadEtVBFCut=15.;
      subleadEtVHadCut=15.;

      diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
      diphotonVHad_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHadCut, subleadEtVHadCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 

      TLorentzVector* jet1=0;
      TLorentzVector* jet2=0;
      TLorentzVector* jet3=0;


      if(diphotonVBF_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut, jet3);

	//if(jet3)
	  //  cout<<"AAA MARCOMM Outside "<<jet3->Pt()<<endl;

        bool VBFpresel = (highestPtJets.first>=0)&&(highestPtJets.second>=0);

	//taken out from if  
	TLorentzVector diphoton = lead_p4+sublead_p4;
	myAll_Mgg =diphoton.M();
	myAllPtHiggs =diphoton.Pt();
        myVBF_Mgg =diphoton.M();

	myAll_phoet1=lead_p4.Et();
	myAll_phoet2=sublead_p4.Et();
	myAll_phoetom1=lead_p4.Et()/diphoton.M();
	myAll_phoetom2=sublead_p4.Et()/diphoton.M();

	if(myAll_phoet1>55.&&myAll_phoet2>25.) {
	  if(myVBF_Mgg>90.&&myVBF_Mgg<190.) {
	    l.FillHist("run",1, l.run, 1.);
	  }
	}

	myAll_nvtx=l.vtx_std_n;

	myAllLeadJPt = 0.;
	myAllSubJPt = 0.;
	myAllLeadJEta = 0.;
	myAllSubJEta = 0.;
	myAll_Mjj = 0.;
	myAlldEta = 0.;
	myAllZep  = 0.;
	myAlldPhi = 0.;
	
	myVBFLeadJPt = 0.;
	myVBFSubJPt = 0.;
	myVBF_Mjj = 0.;
	myVBFdEta = 0.;
	myVBFZep  = 0.;
	myVBFdPhi = 0.;

	if(highestPtJets.first>=0) {

	  jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);

          myAllLeadJPt = jet1->Pt();
          myAllLeadJEta = jet1->Eta();

          myVBFLeadJPt = jet1->Pt();
	}

        if(VBFpresel){

	  //cout<<"MARCOMM "<<highestPtJets.first<<" "<<highestPtJets.second<<endl;

          jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
          
          //myAllLeadJPt = jet1->Pt();
          myAllSubJPt = jet2->Pt();
          //myAllLeadJEta = jet1->Eta();
          myAllSubJEta = jet2->Eta();
          myAll_Mjj = dijet.M();
          myAlldEta = fabs(jet1->Eta() - jet2->Eta());
          myAllZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myAlldPhi = fabs(diphoton.DeltaPhi(dijet));

          //myVBFLeadJPt = jet1->Pt();
          myVBFSubJPt = jet2->Pt();
          myVBF_Mjj = dijet.M();
          myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
          myVBFZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVBFdPhi = fabs(diphoton.DeltaPhi(dijet));
	}

	//should be done for both earlier
	diphoton_index = std::make_pair( l.dipho_leadind[diphotonVBF_id],  l.dipho_subleadind[diphotonVBF_id] );
	float evweight = newweight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight * pileupWeight;
	if(fabs((float) newweight*pileupWeight-weight)/((float) newweight*pileupWeight+weight)>0.0001) cout<<"################ "<<newweight*pileupWeight<<" "<<weight<<" "<<newweight<<" "<<pileupWeight<<endl;
	//if(newweight*pileupWeight != weight) cout<<"AAA################ "<<newweight*pileupWeight<<" "<<weight<<" "<<newweight<<" "<<pileupWeight<<endl;
	float myweight=1.;
	if(evweight*newweight!=0) myweight=evweight/newweight;

	evweightvbf=evweight;
	myweightvbf=myweight;


	//cout<<" Weights: weight "<<weight<<" newtimepileup" <<newweight*pileupWeight<<" genwei "<<genLevWeight<<" PTHihhs "<<myAllPtHiggs<<""<<genLevWeight<<endl;
	
	int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,0.,2,2,1);
	int category2 = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,0.,2,1,1);

	if(!l.GetCutValue("optree")) {
	  //FILL OPTREE FOR ME

	  int itypepass=cur_type;
	


	  //if(jet3)
	  //cout<<"AAA MARCOMM bef set "<<jet3->Pt()<<endl;

	  SetOutputNtupleVariables(jentry, itypepass, l.dipho_leadind[diphotonVBF_id], l.dipho_subleadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], diphoton.M(), &lead_p4, &sublead_p4, evweight, pileupWeight, jet1, jet2, jet3);
	  
	  optree->Fill();
	  int noptree=optree->GetEntries();
	  if(noptree<1000) cout<<" filling optree n="<<noptree<<endl;
	  if(noptree%1000==0||noptree%1000==1) cout<<" optree n="<<noptree<<endl;
	}

	if(l.ApplyCut("massrange",myAll_Mgg,0)) {

	  l.ApplyCutsFill(0,3,evweight, myweight);
	  //a
	  VBFevent = l.ApplyCutsFill(0,1,evweight, myweight);

	  if(passincl&&!VBFevent) {
	    l.ApplyCutsFill(category,45,evweight, myweight);
	  }

	  //VBFevent = l.ApplyCutsFill(0,1,evweight, evweight);
	  l.ApplyCutsFill(0,5,evweight, myweight);
	  l.ApplyCutsFill(0,30,evweight, myweight);
	  l.ApplyCutsFill(0,31,evweight, myweight);

	  int catrun=0;

	  //vtx_std_n

	  if(myAll_nvtx>9) catrun=1;
	  l.ApplyCutsFill(catrun,11,evweight, myweight);

	  /*
	  if(cur_type==0) {
	    if(l.run>175830) catrun=1;

	    l.ApplyCutsFill(catrun,11,evweight, myweight);
	  }
	  else {
	    l.ApplyCutsFill(0,11,evweight*2.247/4.781,myweight*2.247/4.781);
	    l.ApplyCutsFill(1,11,evweight*2.534/4.781,myweight*2.534/4.781);
	  }
	  */


	  
	  if(VBFevent&&myAll_Mgg>100&&myAll_Mgg<180) {
	    l.ApplyCutsFill(category,44,evweight, myweight);
	    l.ApplyCutsFill(category2,42,evweight, myweight);
	  }
	}



	ccat2=category;
	//if(ccat1!=ccat2) cout<<" MARCO DIFF NUMB OF CAT "<<ccat1<<" "<<ccat2<<endl;

	//cout<<cur_type<<" Final events r,l,e"<< l.run << " " << l.lumis << " " << l.event <<" "<<category2<<" "<<diphoton.M()<<" "<<weight<<endl;

	if (cur_type==0){
	  int faketype=999;
	  eventListText << setprecision(4) <<"Type = "<< faketype <<  " Run = " << l.run << "  LS = " << l.lumis << "  Event = " << l.event << "  SelVtx = " << l.vtx_std_sel << "  CAT4 = " << category % 4 << "  ggM = " << diphoton.M() << " gg_Pt =  " << diphoton.Pt();
	  eventListText << endl;
	}


	if(VBFevent){
	  if (cur_type==0){

          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);


	  if(ccat1!=ccat2) cout<<" MARCO NEW DIFF NUMB OF CAT "<<ccat1<<" "<<ccat2<<endl;

	    eventListTextVBF << setprecision(4) <<"Type = "<< cur_type <<  " Run = " << l.run << "  LS = " << l.lumis << "  Event = " << l.event << "  SelVtx = " << l.vtx_std_sel << "  CAT4 = " << category % 4 << "  ggM = " << diphoton.M() << " gg_Pt =  " << diphoton.Pt();

	    eventListTextVBF << setprecision(4) <<  "Run = " << l.run << "  LS = " << l.lumis <<
              "  Event = " << l.event << "  SelVtx = " << l.dipho_vtxind[diphotonVHad_id]
			     << "  CAT4 = " << category % 4 << "  ggM = " << myVHad_Mgg << " ggPt =  " << diphoton.Pt()
			     << "  jetEta1 = " << jet1->Eta() << "  jetEta2 = " << jet2->Eta() 
			     << "  jetPhi1 = " << jet1->Phi() << "  jetPhi2 = " << jet2->Phi()
			     <<  "  jetEt1 = " << jet1->Et() << "  jetEt2 = "  << jet2->Et()
			     << " Mjj " << myAll_Mjj
			     << " dEtajj " << myAlldEta
			     << " Zeppenfeld " << myAllZep
			     << " dPhijjgg " << myAlldPhi << " VH itype " <<cur_type << endl;
	      
	    eventListTextVBF << setprecision(4) <<"Type = "<< cur_type <<  "Run = " << l.run << "  LS = " << l.lumis << "  Event = " << l.event << "  SetchangedVtx = " << l.vtx_std_sel << "  CAT4 = " << category % 4 << "  ggM = " << diphoton.M() << " gg_Pt =  " << diphoton.Pt();
	    eventListTextVBF << endl;
	    eventListTextVBF << setprecision(4) <<" phoet   "<<myAll_phoet1<<" "<<myAll_phoet2<<endl;
	    eventListTextVBF << setprecision(4) <<" phoetom "<<myAll_phoetom1<<" "<<myAll_phoetom2<<endl;
	    eventListTextVBF << setprecision(4) <<" jetet   "<<myAllLeadJPt<<" "<<myAllSubJPt<<endl;
	    eventListTextVBF << setprecision(4) <<" jeteta  "<<myAllLeadJEta<<" "<<myAllSubJEta<<endl;
	    eventListTextVBF << setprecision(4) <<" jets deta "<<myAlldEta<<" zep "<<myAllZep<<" M "<<myAll_Mjj<<" dphi "<<myAlldPhi <<endl;
	    eventListTextVBF << setprecision(4) <<" phoeta   "<<lead_p4.Eta()<<" "<<sublead_p4.Eta()<<endl;
	    eventListTextVBF << setprecision(4) <<" phophi   "<<lead_p4.Phi()<<" "<<sublead_p4.Phi()<<endl;
			    
	  }
	}
	

      	if(l.ApplyCut("massrange",myAll_Mgg,0)) {

	  VHadevent = 0; //l.ApplyCutsFill(0,2,evweight, myweight);

	//l.ApplyCutsFill(0,6,evweight, myweight);
	}

      }
      /*
      if(diphotonVHad_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );

        bool VHadpresel = (highestPtJets.first>=0)&&(highestPtJets.second>=0);
  
	//taken out from if  
	TLorentzVector diphoton = lead_p4+sublead_p4;
        myVHad_Mgg =diphoton.M();
	
	myVHadLeadJPt = 0.;
	myVHadSubJPt = 0.;
	myVHad_Mjj = 0.;
	myVHaddEta = 0.;
	myVHadZep  = 0.;
	myVHaddPhi = 0.;

	if(highestPtJets.first>=0) {

	  TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);

          myVHadLeadJPt = jet1->Pt();
	}

	if(VHadpresel){
          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
	  
          TLorentzVector diphoton = lead_p4+sublead_p4;
          
          //myVHadLeadJPt = jet1->Pt();
          myVHadSubJPt = jet2->Pt();
          myVHad_Mjj = dijet.M();
          myVHaddEta = fabs(jet1->Eta() - jet2->Eta());
          myVHadZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVHaddPhi = fabs(diphoton.DeltaPhi(dijet));
	}

	float evweight = newweight * smeared_pho_weight[l.dipho_leadind[diphotonVHad_id]] * smeared_pho_weight[l.dipho_subleadind[diphotonVHad_id]] * genLevWeight * pileupWeight;
	if(fabs((float) newweight*pileupWeight-weight)/((float) newweight*pileupWeight+weight)>0.0001) cout<<"################ "<<newweight*pileupWeight<<" "<<weight<<" "<<newweight<<" "<<pileupWeight<<endl;
	//if(newweight*pileupWeight != weight) cout<<"BBB################ "<<newweight*pileupWeight<<" "<<weight<<" "<<newweight<<" "<<pileupWeight<<endl;
	float myweight=1.;
	if(evweight*newweight!=0) myweight=evweight/newweight;
	
	VHadevent = l.ApplyCutsFill(0,2,evweight, myweight);
	l.ApplyCutsFill(0,6,evweight, myweight);
       
      }
      */
    }



    if(VBFevent) {
      if(l.ApplyCut("massrange",myAll_Mgg,0)) {
	l.FillHist("Mass5cat",0,myAll_Mgg,evweightvbf);
      }
      l.FillCounter("Mass5cat",myweightvbf,0);
    }
    else if(passincl) {
      if(l.ApplyCut("massrange",myAll_Mgg,0)) {
	l.FillHist("Mass5cat",catincl+1,myAll_Mgg,evweightincl);
      }
      l.FillCounter("Mass5cat",myweightincl,catincl+1);
    }
    
    if(includeVBF&&VBFevent) diphoton_id = diphotonVBF_id;
    else if(includeVHad&&VHadevent) diphoton_id = diphotonVHad_id;
    
    // CP

    /// std::cerr << "Selected pair " << l.dipho_n << " " << diphoton_id << std::endl;
    if (diphoton_id > -1 ) {

	diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
    	// bring all the weights together: lumi & Xsection, single gammas, pt kfactor
	float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

	l.countersred[diPhoCounter_]++;

	TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
	float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
	TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

	bool CorrectVertex;
	
  
  // FIXME pass smeared R9
	int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
	int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
	if( cur_type != 0 && doMCSmearing ) {
	    float pth = Higgs.Pt();
	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
		float rewei=1.;
		(*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
		if( rewei < 0. ) {
		    std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		    assert(0);
		}
		evweight *= rewei;
	    }
	    CorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
	}
	float mass    = Higgs.M();
	float ptHiggs = Higgs.Pt();
      
	assert( evweight >= 0. ); 

	l.FillCounter( "Accepted", weight );
	l.FillCounter( "Smeared", evweight );
	sumaccept += weight;
 	sumsmear += evweight;

	if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
	else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;



	// control plots 
	l.FillHist("all_mass",0, Higgs.M(), evweight);
	l.FillHist("all_mass",category+1, Higgs.M(), evweight);
	if( mass>=massMin && mass<=massMax  ) {
		l.FillHist("mass",0, Higgs.M(), evweight);
		l.FillHist("pt",0, Higgs.Pt(), evweight);
		l.FillHist("eta",0, Higgs.Eta(), evweight);
		
		l.FillHist("pho_pt",0,lead_p4.Pt(), evweight);
		l.FillHist("pho1_pt",0,lead_p4.Pt(), evweight);
		l.FillHist("pho_eta",0,lead_p4.Eta(), evweight);
		l.FillHist("pho1_eta",0,lead_p4.Eta(), evweight);
		l.FillHist("pho_r9",0, lead_r9, evweight);
		l.FillHist("pho1_r9",0, lead_r9, evweight);
		
		l.FillHist("pho_pt",0,sublead_p4.Pt(), evweight);
		l.FillHist("pho2_pt",0,sublead_p4.Pt(), evweight);
		l.FillHist("pho_eta",0,sublead_p4.Eta(), evweight);
		l.FillHist("pho2_eta",0,sublead_p4.Eta(), evweight);
		l.FillHist("pho_r9",0, sublead_r9, evweight);
		l.FillHist("pho1_r9",0, sublead_r9, evweight);
		
		l.FillHist("mass",category+1, Higgs.M(), evweight);
		l.FillHist("pt",category+1, Higgs.Pt(), evweight);
		l.FillHist("eta",category+1, Higgs.Eta(), evweight);
		
		l.FillHist("pho_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho1_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,lead_p4.Eta(), evweight);
		l.FillHist("pho1_eta",category+1,lead_p4.Eta(), evweight);
		l.FillHist("pho_r9",category+1, lead_r9, evweight);
		l.FillHist("pho1_r9",category+1, lead_r9, evweight);
		
		l.FillHist("pho_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho2_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,sublead_p4.Eta(), evweight);
		l.FillHist("pho2_eta",category+1,sublead_p4.Eta(), evweight);
		l.FillHist("pho_r9",category+1, sublead_r9, evweight);
		l.FillHist("pho1_r9",category+1, sublead_r9, evweight);
		
		l.FillHist("pho_n",category+1,l.pho_n, evweight);
	}

	if (cur_type==0){
	  eventListText << setprecision(4) <<"Type = "<< cur_type <<  " Run = " << l.run << "  LS = " << l.lumis << "  Event = " << l.event << "  SelVtx = " << l.vtx_std_sel << "  CAT4 = " << category % 4 << "  ggM = " << mass << " gg_Pt =  " << ptHiggs;
	  eventListText << endl;
	}
       
	// --------------------------------------------------------------------------------------------- 
	if (cur_type == 0 ){
	    l.rooContainer->InputDataPoint("data_mass",category,mass);
	}
	if (cur_type > 0 && cur_type != 3 && cur_type != 4)
	    l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
	else if (cur_type < 0){
	    l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type),category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_rv",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type)+"_wv",category,mass,evweight);
	}
       
    }
   
   
    // Systematics
    if( cur_type != 0 && doMCSmearing ) { 
	// fill steps for syst uncertainty study
	float systStep = systRange / (float)nSystSteps;
	// di-photon smearers systematics
	if (diphoton_id > -1 ) {
	       
	    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
	 
	    for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
		std::vector<double> mass_errors;
		std::vector<double> weights;
		std::vector<int>    categories;
	   
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	     
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
		    double genLevWeightSyst=1; 
	     
		    for(std::vector<BaseGenLevelSmearer *>::iterator sj=genLevelSmearers_.begin(); sj!= genLevelSmearers_.end(); ++sj ) {
			float swei=1.;
			if( *si == *sj ) { 
			    (*si)->smearEvent(swei, gP4, l.pu_n, cur_type, syst_shift );
			} else {
			    (*sj)->smearEvent(swei, gP4, l.pu_n, cur_type, 0. );
			}
			genLevWeightSyst *= swei;
		    }
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeightSyst;
	     
		    float mass = Higgs.M();
        
        if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
        else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;
		    categories.push_back(category);
		    mass_errors.push_back(mass);
		    weights.push_back(evweight);
		}// end loop on systematics steps
	   
		if (cur_type < 0){
		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
		}
	    }// end loop on smearers 
		 

	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		std::vector<double> mass_errors;
		std::vector<double> weights;
		std::vector<int> categories;
		       
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
			       
		    // restart with 'fresh' wait for this round of systematics
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
			       
		    // FIXME pass smeared R9 and di-photon
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
		    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		    for(std::vector<BaseDiPhotonSmearer *>::iterator sj=diPhotonSmearers_.begin(); sj!= diPhotonSmearers_.end(); ++sj ) {
			float swei=1.;
			float pth = Higgs.Pt();
			if( *si == *sj ) { 
			    (*si)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), syst_shift );
			} else { 
			    (*sj)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			}
			evweight *= swei;
		    }
		    float mass = Higgs.M();
        if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
        else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;
		    categories.push_back(category);
		    mass_errors.push_back(mass);
		    weights.push_back(evweight);
		}
		if (cur_type < 0){
		    l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
		}
	    }

	}
       
	// loop over the smearers included in the systematics study
	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	    std::vector<double> mass_errors;
	    std::vector<double> weights;
	    std::vector<int> categories;
	   
	    // loop over syst shift
	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		if( syst_shift == 0. ) { continue; } // skip the central value
		// smear the photons 
		for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		    std::vector<std::vector<bool> > p;
		    //std::cout << "GF check: " <<  l.pho_residCorrEnergy[ipho] << "  " << l.pho_residCorrResn[ipho] << std::endl;
		    PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
						/// *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
						((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
						energyCorrected[ipho],
						l.pho_isEB[ipho], l.pho_r9[ipho],
						l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_));
		   
		    float pweight = 1.;
		    for(std::vector<BaseSmearer *>::iterator  sj=photonSmearers_.begin(); sj!= photonSmearers_.end(); ++sj ) {
			float sweight = 1.;
			if( *si == *sj ) {
			    // move the smearer under study by syst_shift
			    (*si)->smearPhoton(phoInfo,sweight,l.run,syst_shift);
			} else {
			    // for the other use the nominal points
			    (*sj)->smearPhoton(phoInfo,sweight,l.run,0.);
			}
			pweight *= sweight;
		    }
		    smeared_pho_energy[ipho] = phoInfo.energy();
		    smeared_pho_r9[ipho] = phoInfo.r9();
		    smeared_pho_weight[ipho] = pweight;
		}
	       
		// analyze the event
		// FIXME pass smeared R9
		int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
    int diphotonVBF_id = -1;
    int diphotonVHad_id = -1;
    bool VBFevent = false;
    bool VHadevent = false;
    std::pair<int,int> highestPtJets(-1,-1);

    if((includeVBF || includeVHad)&&l.jet_algoPF1_n>1) {

      int applyPtoverM=0;

      if(includeVBF) diphotonVBF_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
      if(includeVHad) diphotonVHad_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHadCut, subleadEtVHadCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 

      if(diphotonVBF_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
	      TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVBF_id], l.dipho_vtxind[diphotonVBF_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );

        bool VBFpresel = (highestPtJets.first>=0)&&(highestPtJets.second>=0);
  

        if(VBFpresel){
          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
	        
          TLorentzVector diphoton = lead_p4+sublead_p4;
          
          myVBFLeadJPt = jet1->Pt();
          myVBFSubJPt = jet2->Pt();
          myVBF_Mjj = dijet.M();
          myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
          myVBFZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVBFdPhi = fabs(diphoton.DeltaPhi(dijet));
          myVBF_Mgg =diphoton.M();

          VBFevent = l.ApplyCuts(0,1);
        }
      }
      
      if(diphotonVHad_id>-1){
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
	      TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHad_id], l.dipho_vtxind[diphotonVHad_id], &smeared_pho_energy[0]);
        float jet1ptcut =0.0;
        float jet2ptcut =0.0;
        
        highestPtJets = Select2HighestPtJets(l, lead_p4, sublead_p4, jet1ptcut, jet2ptcut );

        bool VHadpresel = (highestPtJets.first>=0)&&(highestPtJets.second>=0);
  
        if(VHadpresel){
          TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
          TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
          TLorentzVector dijet = (*jet1) + (*jet2);
	        
          TLorentzVector diphoton = lead_p4+sublead_p4;
          
          myVHadLeadJPt = jet1->Pt();
          myVHadSubJPt = jet2->Pt();
          myVHad_Mjj = dijet.M();
          myVHaddEta = fabs(jet1->Eta() - jet2->Eta());
          myVHadZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
          myVHaddPhi = fabs(diphoton.DeltaPhi(dijet));
          myVHad_Mgg =diphoton.M();

          VHadevent = 0; //l.ApplyCuts(0,2);
        }
      }
    }
    
    if(VBFevent) diphoton_id = diphotonVBF_id;
    else if(VHadevent) diphoton_id = diphotonVHad_id;
	       
		if (diphoton_id > -1 ) {
		   
		    diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] *genLevWeight;
		   
		    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
		   
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
		    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		    if( cur_type != 0 && doMCSmearing ) {
			for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
			    float rewei=1.;
			    float pth = Higgs.Pt();
			    (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			    evweight *= rewei;
			}
		    }
		    float mass = Higgs.M();
		   
              if(VBFevent) category=nEtaCategories*nR9Categories*nPtCategories;
              else if(VHadevent) category=nEtaCategories*nR9Categories*nPtCategories+1;
	       	    categories.push_back(category);
	            mass_errors.push_back(mass);
	            weights.push_back(evweight);

		} else {
		    mass_errors.push_back(0.);   
		    weights.push_back(0.);   
		    categories.push_back(-1);
		}
	       
	    }
	    if (cur_type < 0){
	       l.rooContainer->InputSystematicSet("sig_"+GetSignalLabel(cur_type),(*si)->name(),categories,mass_errors,weights);
	    }
       
	}
       
       
    }
   
    if(PADEBUG) 
	cout<<"myFillHistRed END"<<endl;
}

// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
void StatAnalysisExclusive::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool StatAnalysisExclusive::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}
// ----------------------------------------------------------------------------------------------------
double StatAnalysisExclusive::GetDifferentialKfactor(double gPT, int Mass)
{

/*  
    if (Mass <=110 ) return thm110->GetBinContent(thm110->FindFixBin(gPT));
    else if (Mass ==120 ) return thm120->GetBinContent(thm120->FindFixBin(gPT));
    else if (Mass ==130 ) return thm130->GetBinContent(thm130->FindFixBin(gPT));
    else if (Mass ==140 ) return thm140->GetBinContent(thm140->FindFixBin(gPT));
    else if (Mass ==115 ) return (0.5*thm110->GetBinContent(thm110->FindFixBin(gPT)) +0.5*thm120->GetBinContent(thm120->FindFixBin(gPT)));
*/
    return 1.0;
/*
  int  genMasses[4] = {110,120,130,140};
  if (Mass<=genMasses[0] ) return kfactorHistograms[0]->GetBinContent(kfactorHistograms[0]->FindBin(gPT));
  else if (Mass<genMasses[nMasses-1]) {

  TH1D *hm1,*hm2;
  double m1=0,m2=0;
  for (int m=0;m<nMasses;m++){
  if (Mass<genMasses[m+1]){
  hm1=kfactorHistograms[m];
  hm2=kfactorHistograms[m+1];
  m1 = genMasses[m];
  m2 = genMasses[m+1];
  //	cout << "Gen Mass: "<< Mass << " Using "<<m1<< " " << m2<< " Hist name check " << hm1->GetName()<<" " <<hm2->GetName()<<endl;
  break;
  }
  }
  if ((int)Mass == (int)m1 ){
  //cout << "Found the appropriate historgam "<<hm1->GetName()<<endl;
  return hm1->GetBinContent(hm1->FindBin(gPT));
  } else {

  TH1D *hm = (TH1D*) hm1->Clone("hm");
  double alpha = ((float) (Mass-m1))/(m2-m1); // make sure ms are not integers
  hm->Add(hm1,hm2,alpha,(1-alpha));
  return hm->GetBinContent(hm->GetBinContent(hm->FindBin(gPT)));
  }

  }
  else return kfactorHistograms[nMasses-1]->GetBinContent(kfactorHistograms[nMasses-1]->FindBin(gPT));
*/
}

void StatAnalysisExclusive::FillSignalLabelMap(){

	// Basically A Map of the ID (type) to the signal's name which can be filled Now:
  signalLabels[-57]="ggh_mass_m123";
  signalLabels[-58]="vbf_mass_m123";
  signalLabels[-60]="wzh_mass_m123";
  signalLabels[-59]="tth_mass_m123";
  signalLabels[-53]="ggh_mass_m121";
  signalLabels[-54]="vbf_mass_m121";
  signalLabels[-56]="wzh_mass_m121";
  signalLabels[-55]="tth_mass_m121";
  signalLabels[-65]="ggh_mass_m160";
  signalLabels[-66]="vbf_mass_m160";
  signalLabels[-68]="wzh_mass_m160";
  signalLabels[-67]="tth_mass_m160";
  signalLabels[-61]="ggh_mass_m155";
  signalLabels[-62]="vbf_mass_m155";
  signalLabels[-64]="wzh_mass_m155";
  signalLabels[-63]="tth_mass_m155";
  signalLabels[-49]="ggh_mass_m150";
  signalLabels[-50]="vbf_mass_m150";
  signalLabels[-52]="wzh_mass_m150";
  signalLabels[-51]="tth_mass_m150";
  signalLabels[-45]="ggh_mass_m145";
  signalLabels[-46]="vbf_mass_m145";
  signalLabels[-48]="wzh_mass_m145";
  signalLabels[-47]="tth_mass_m145";
  signalLabels[-33]="ggh_mass_m140";
  signalLabels[-34]="vbf_mass_m140";
  signalLabels[-36]="wzh_mass_m140";
  signalLabels[-35]="tth_mass_m140";
  signalLabels[-41]="ggh_mass_m135";
  signalLabels[-42]="vbf_mass_m135";
  signalLabels[-44]="wzh_mass_m135";
  signalLabels[-43]="tth_mass_m135";
  signalLabels[-29]="ggh_mass_m130";
  signalLabels[-30]="vbf_mass_m130";
  signalLabels[-32]="wzh_mass_m130";
  signalLabels[-31]="tth_mass_m130";
  signalLabels[-37]="ggh_mass_m125";
  signalLabels[-38]="vbf_mass_m125";
  signalLabels[-40]="wzh_mass_m125";
  signalLabels[-39]="tth_mass_m125";
  signalLabels[-25]="ggh_mass_m120";
  signalLabels[-26]="vbf_mass_m120";
  signalLabels[-28]="wzh_mass_m120";
  signalLabels[-27]="tth_mass_m120";
  signalLabels[-21]="ggh_mass_m115";
  signalLabels[-22]="vbf_mass_m115";
  signalLabels[-24]="wzh_mass_m115";
  signalLabels[-23]="tth_mass_m115";
  signalLabels[-17]="ggh_mass_m110";
  signalLabels[-18]="vbf_mass_m110";
  signalLabels[-19]="wzh_mass_m110";
  signalLabels[-20]="tth_mass_m110";
  signalLabels[-13]="ggh_mass_m105";
  signalLabels[-14]="vbf_mass_m105";
  signalLabels[-16]="wzh_mass_m105";
  signalLabels[-15]="tth_mass_m105";
  signalLabels[-69]="ggh_mass_m100";
  signalLabels[-70]="vbf_mass_m100";
  signalLabels[-72]="wzh_mass_m100";
  signalLabels[-71]="tth_mass_m100";
}

std::string StatAnalysisExclusive::GetSignalLabel(int id){
	
	// For the lazy man, can return a memeber of the map rather than doing it yourself
	std::map<int,std::string>::iterator it = signalLabels.find(id);

	if (it!=signalLabels.end()){
		return it->second;
		
	} else { 

		std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
		return "NULL";
	}
	
}

std::pair<int, int> StatAnalysisExclusive::Select2HighestPtJets(LoopAll& l, TLorentzVector& leadpho, TLorentzVector& subleadpho, float jtLMinPt, float jtTMinPt, TLorentzVector* jet3){

  std::pair<int, int> myJets(-1,-1);
  std::pair<int, int> fail(-1,-1);

  std::pair<int, int> myJetsnew(-1,-1);
  std::pair<float, float> myJetspt(-1.,-1.);
  std::pair<float, float> myJetsnewpt(-1.,-1.);

  float dr2pho = 0.5;
  float dr2jet = 0.5; //Shoud change marco

  TLorentzVector* j1p4;
  TLorentzVector* j2p4;
  float j1pt=-1;
  float j2pt=-1;

  /*
  // select ighest pt jets
  // veto jets close to photons or each other
  for(int j1_i=0; j1_i<l.jet_algoPF1_n; j1_i++){
    j1p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(j1_i);
    if(fabs(j1p4->Eta()) > 4.7) continue;
    if(j1p4->DeltaR(leadpho) < dr2pho) continue;
    if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
    j1pt=j1p4->Pt();
    if(j1pt<jtTMinPt) continue; //jtT
    for(int j2_i=j1_i+1; j2_i<l.jet_algoPF1_n; j2_i++){
      j2p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(j2_i);
      if(fabs(j2p4->Eta()) > 4.7) continue;
      if(j2p4->DeltaR(leadpho) < dr2pho) continue;
      if(j2p4->DeltaR(subleadpho) < dr2pho) continue;
      if(j2p4->DeltaR(*j1p4) < dr2jet) continue;
      j2pt=j2p4->Pt();
      
      if(j2pt<jtTMinPt) continue; //jtT
      if(std::max(j1pt,j2pt)<jtLMinPt) continue;

      if(j1pt>j2pt){
        jtLMinPt=j1pt; //??? why using something for something else???
        jtTMinPt=j2pt; //??? why using something for something else???

        myJets.first = j1_i;
        myJets.second = j2_i;
      } else {
        jtLMinPt=j2pt;
        jtTMinPt=j1pt;

        myJets.first = j2_i;
        myJets.second = j1_i;
      }
    }
  }
  */

  float pt3=0;
  int ind3=-1;

  for(int j1_i=0; j1_i<l.jet_algoPF1_n; j1_i++){
    j1p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(j1_i);
    if(fabs(j1p4->Eta()) > 4.7) continue;
    if(j1p4->DeltaR(leadpho) < dr2pho) continue;
    if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
    j1pt=j1p4->Pt();

    //cout<<"AAA MARCOMM "<<j1_i<<" "<<j1p4->Pt()<<" "<<j1p4->Eta()<<endl;

    //if(j1pt<jtTMinPt) continue;

    if(j1pt>myJetsnewpt.first) {

      pt3=myJetsnewpt.second;
      ind3=myJetsnew.second;

      myJetsnew.second=myJetsnew.first;
      myJetsnewpt.second=myJetsnewpt.first;
      myJetsnewpt.first=j1pt;
      myJetsnew.first=j1_i;
    }
    else if(j1pt>myJetsnewpt.second) {

      pt3=myJetsnewpt.second;
      ind3=myJetsnew.second;

      myJetsnewpt.second=j1pt;
      myJetsnew.second=j1_i;
    }
    else if(j1pt>pt3) {

      pt3=j1pt;
      ind3=j1_i;

    }
  }

  if(ind3>-1) 
    jet3 = (TLorentzVector*) l.jet_algoPF1_p4->At(ind3);

  //cout<<"AAA MARCOMM "<<l.jet_algoPF1_n<<" "<<myJetsnew.first<<" "<<myJetsnew.second<<endl;

  //if(jet3)
  //cout<<"AAA MARCOMM "<<ind3<<" "<<jet3->Pt()<<endl;

  //if(myJets.first==-1) return fail;
  //return myJets;

  /*
  if(myJetsnew.second!=-1&&myJetsnewpt.first>jtTMinPt&&myJetsnewpt.second>jtTMinPt) {
    if(myJetsnew!=myJets) {
      j1p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(myJetsnew.first);
      j2p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(myJetsnew.second);
      float dr=j2p4->DeltaR(*j1p4);
      //cout<<"myJetsnew myJets "<<myJetsnew.first<<myJetsnew.second<<myJets.first<<myJets.second<<" dr "<<dr<<endl;
      //cout<<"myJetsnew myJets "<<myJetspt.first<<" "<<myJetspt.second<<" "<<jtLMinPt<<jtTMinPt<<endl;
    }
  }
  */

  return myJetsnew;
}



int  StatAnalysisExclusive::RescaleJetEnergy(LoopAll& l) {
  for (int i = 0; i<l.jet_algoPF1_n; i++) {
    TLorentzVector * thisjet = (TLorentzVector *) l.jet_algoPF1_p4->At(i);
    *thisjet*=l.jet_algoPF1_erescale[i];
  }
  return 1;
}

void StatAnalysisExclusive::HggBookOptree() {

  //opfile = new TFile("optNtuple.root","RECREATE","optimization ntuple");
  //opfile->cd();
  //ll->outputFile->cd();  //MARCO CHECK THIS
  ll->hfilereal->cd();  //MARCO CHECK THIS
  optree= new TTree("ntuple","Hgg optimization Tree");

  optree->Branch("j1pt",&t_j1pt,"j1pt/F",2000000);
  optree->Branch("j1eta",&t_j1eta,"j1eta/F",2000000);
  optree->Branch("j1phi",&t_j1phi,"j1phi/F",2000000);
  optree->Branch("j2pt",&t_j2pt,"j2pt/F",2000000);
  optree->Branch("j2eta",&t_j2eta,"j2eta/F",2000000);
  optree->Branch("j2phi",&t_j2phi,"j2phi/F",2000000);

  optree->Branch("j3pt",&t_j3pt,"j3pt/F",2000000);
  optree->Branch("j3eta",&t_j3eta,"j3eta/F",2000000);
  optree->Branch("j3phi",&t_j3phi,"j3phi/F",2000000);

  optree->Branch("jjmass",&t_jjmass,"jjmass/F",2000000);
  optree->Branch("jjzep",&t_jjzep,"jjzep/F",2000000);
  optree->Branch("jjdeta",&t_jjdeta,"jjdeta/F",2000000);
  optree->Branch("jjdphi",&t_jjdphi,"jjdphi/F",2000000);


  //event vars
  optree->Branch("run",&t_run,"run/I",2000000);
  optree->Branch("lumis",&t_lumis,"lumis/I",2000000);
  optree->Branch("event",&t_event,"event/I",2000000);
  optree->Branch("itype",&t_itype,"itype/I",2000000);
  optree->Branch("processid",&t_processid,"processid/I",2000000);
  optree->Branch("w",&t_w,"w/F",2000000);
  optree->Branch("wpu",&t_wpu,"wpu/F",2000000);
  //diphoton variables
  optree->Branch("mass",&t_mass,"mass/F",2000000);
  optree->Branch("dmom0",&t_dmom0,"dmom0/F",2000000);
  optree->Branch("dmom",&t_dmom,"dmom/F",2000000);
  optree->Branch("deltaM",&t_deltaM,"deltaM/F",2000000);
  optree->Branch("category",&t_category,"category/I",2000000);
  optree->Branch("diphocat2r92eta",&t_diphocat2r92eta,"diphocat2r92eta/I",2000000);
  optree->Branch("diphosubcat4",&t_diphosubcat4,"diphosubcat4/I",2000000);
  optree->Branch("barrel",&t_barrel,"barrel/I",2000000);
  optree->Branch("diphor9",&t_diphor9,"diphor9/F",2000000);
  optree->Branch("diphoeta",&t_diphoeta,"diphoeta/F",2000000);
  optree->Branch("costhetastar",&t_costhetastar,"costhetastar/F",2000000);
  optree->Branch("diphopt",&t_diphopt,"diphopt/F",2000000);
  optree->Branch("diphopz",&t_diphopz,"diphopz/F",2000000);
  optree->Branch("deltar",&t_deltar,"deltar/F",2000000);
  optree->Branch("etamax",&t_etamax,"etamax/F",2000000);
  optree->Branch("etamin",&t_etamin,"etamin/F",2000000);
  optree->Branch("deta",&t_deta,"deta/F",2000000);
  //lead photon variables
  optree->Branch("leadcat",&t_leadcat,"leadcat/I",2000000);
  optree->Branch("leadr9",&t_leadr9,"leadr9/F",2000000);
  optree->Branch("leadeta",&t_leadeta,"leadeta/F",2000000);
  optree->Branch("leadpt",&t_leadpt,"leadpt/F",2000000);
  optree->Branch("leadgenmatch",&t_leadgenmatch,"leadgenmatch/I",2000000);
  optree->Branch("leadbarrel",&t_leadbarrel,"leadbarrel/I",2000000);
  optree->Branch("leadsee",&t_leadsee,"leadsee/F",2000000);
  optree->Branch("leadpi0nn",&t_leadpi0nn,"leadpi0nn/F",2000000);
  //sublead photon variables
  optree->Branch("subleadcat",&t_subleadcat,"subleadcat/I",2000000);
  optree->Branch("subleadr9",&t_subleadr9,"subleadr9/F",2000000);
  optree->Branch("subleadeta",&t_subleadeta,"subleadeta/F",2000000);
  optree->Branch("subleadpt",&t_subleadpt,"subleadpt/F",2000000);
  optree->Branch("subleadgenmatch",&t_subleadgenmatch,"subleadgenmatch/I",2000000);
  optree->Branch("subleadbarrel",&t_subleadbarrel,"subleadbarrel/I",2000000);
  optree->Branch("subleadsee",&t_subleadsee,"subleadsee/F",2000000);
  optree->Branch("subleadpi0nn",&t_subleadpi0nn,"subleadpi0nn/F",2000000);
  //PTDR isolation vars
  optree->Branch("leadtrkiso",&t_leadtrkiso,"leadtrkiso/I",2000000);
  optree->Branch("leadecaliso",&t_leadecaliso,"leadecaliso/F",2000000);
  optree->Branch("leadhcaliso",&t_leadhcaliso,"leadhcaliso/F",2000000);
  optree->Branch("leadhovere",&t_leadhovere,"leadhovere/F",2000000);
  optree->Branch("subleadtrkiso",&t_subleadtrkiso,"subleadtrkiso/I",2000000);
  optree->Branch("subleadecaliso",&t_subleadecaliso,"subleadecaliso/F",2000000);
  optree->Branch("subleadhcaliso",&t_subleadhcaliso,"subleadhcaliso/F",2000000);
  optree->Branch("subleadhovere",&t_subleadhovere,"subleadhovere/F",2000000);
  //lead tracker

  /*       jgb
  optree->Branch("leadtrkdeltar0",&t_leadtrkdeltar0,"leadtrkdeltar0/F",2000000);
  optree->Branch("leadtrkdeltar15",&t_leadtrkdeltar15,"leadtrkdeltar15/F",2000000);
  optree->Branch("leadtrkdeltar20",&t_leadtrkdeltar20,"leadtrkdeltar20/F",2000000);
  optree->Branch("leadtrkdeltar25",&t_leadtrkdeltar25,"leadtrkdeltar25/F",2000000);
  optree->Branch("leadtrkecone20",&t_leadtrkecone20,"leadtrkecone20/F",2000000);
  optree->Branch("leadtrkecone25",&t_leadtrkecone25,"leadtrkecone25/F",2000000);
  optree->Branch("leadtrkecone30",&t_leadtrkecone30,"leadtrkecone30/F",2000000);
  optree->Branch("leadtrkecone35",&t_leadtrkecone35,"leadtrkecone35/F",2000000);
  optree->Branch("leadtrkecone40",&t_leadtrkecone40,"leadtrkecone40/F",2000000);
  optree->Branch("leadtrkecone45",&t_leadtrkecone45,"leadtrkecone45/F",2000000);
  //lead ecal
  optree->Branch("leadecalhitsJ_060330",&t_leadecalhitsJ_060330,"leadecalhitsJ_060330/F",2000000);
  optree->Branch("leadecal_048_32_16_018",&t_leadecal_048_32_16_018,"leadecal_048_32_16_018/F",2000000);
  optree->Branch("leadecal_048_32_14_015",&t_leadecal_048_32_14_015,"leadecal_048_32_14_015/F",2000000);
  optree->Branch("leadecal_048_32_18_018",&t_leadecal_048_32_18_018,"leadecal_048_32_18_018/F",2000000);
  optree->Branch("leadecal_075_36_20_015",&t_leadecal_075_36_20_015,"leadecal_075_36_20_015/F",2000000);
  optree->Branch("leadecal_060_32_20_018",&t_leadecal_060_32_20_018,"leadecal_060_32_20_018/F",2000000);
  optree->Branch("leadecal_075_36_12_015",&t_leadecal_075_36_12_015,"leadecal_075_36_12_015/F",2000000);
  //lead other
  optree->Branch("leadtrkplusecal",&t_leadtrkplusecal,"leadtrkplusecal/F",2000000);
  //sublead tracker
  optree->Branch("subleadtrkdeltar0",&t_subleadtrkdeltar0,"subleadtrkdeltar0/F",2000000);
  optree->Branch("subleadtrkdeltar15",&t_subleadtrkdeltar15,"subleadtrkdeltar15/F",2000000);
  optree->Branch("subleadtrkdeltar20",&t_subleadtrkdeltar20,"subleadtrkdeltar20/F",2000000);
  optree->Branch("subleadtrkdeltar25",&t_subleadtrkdeltar25,"subleadtrkdeltar25/F",2000000);
  optree->Branch("subleadtrkecone20",&t_subleadtrkecone20,"subleadtrkecone20/F",2000000);
  optree->Branch("subleadtrkecone25",&t_subleadtrkecone25,"subleadtrkecone25/F",2000000);
  optree->Branch("subleadtrkecone30",&t_subleadtrkecone30,"subleadtrkecone30/F",2000000);
  optree->Branch("subleadtrkecone35",&t_subleadtrkecone35,"subleadtrkecone35/F",2000000);
  optree->Branch("subleadtrkecone40",&t_subleadtrkecone40,"subleadtrkecone40/F",2000000);
  optree->Branch("subleadtrkecone45",&t_subleadtrkecone45,"subleadtrkecone45/F",2000000);
  //sublead ecal
  optree->Branch("subleadecalhitsJ_060330",&t_subleadecalhitsJ_060330,"subleadecalhitsJ_060330/F",2000000);
  optree->Branch("subleadecal_048_32_16_018",&t_subleadecal_048_32_16_018,"subleadecal_048_32_16_018/F",2000000);
  optree->Branch("subleadecal_048_32_14_015",&t_subleadecal_048_32_14_015,"subleadecal_048_32_14_015/F",2000000);
  optree->Branch("subleadecal_048_32_18_018",&t_subleadecal_048_32_18_018,"subleadecal_048_32_18_018/F",2000000);
  optree->Branch("subleadecal_075_36_20_015",&t_subleadecal_075_36_20_015,"subleadecal_075_36_20_015/F",2000000);
  optree->Branch("subleadecal_060_32_20_018",&t_subleadecal_060_32_20_018,"subleadecal_060_32_20_018/F",2000000);
  optree->Branch("subleadecal_075_36_12_015",&t_subleadecal_075_36_12_015,"subleadecal_075_36_12_015/F",2000000);
  */

  optree->Branch("pvtx", &t_pvtx, "pvtx/F");
  optree->Branch("sigma_mz", &t_sigma_mz, "sigma_mz/F");
  optree->Branch("leadphoidmitmva", &t_leadphoidmitmva, "leadphoidmitmva/F");
  optree->Branch("subleadphoidmitmva", &t_subleadphoidmitmva, "subleadphoidmitmva/F");
  optree->Branch("diphomitmva", &t_diphomitmva, "diphomitmva/F");
  optree->Branch("bsZ", &t_bsZ, "bsZ/F");

  optree->Branch("nvtx",&t_nvtx,"nvtx/F",2000000);
  optree->Branch("rho",&t_rho,"rho/F",2000000);
  optree->Branch("leadcutindex",&t_leadcutindex,"leadcutindex/F",2000000);
  optree->Branch("subleadcutindex",&t_subleadcutindex,"subleadcutindex/F",2000000);
  optree->Branch("leadci6cindex",&t_leadci6cindex,"leadci6cindex/F",2000000);
  optree->Branch("subleadci6cindex",&t_subleadci6cindex,"subleadci6cindex/F",2000000);
  optree->Branch("leadci6cpfindex",&t_leadci6cpfindex,"leadci6cpfindex/F",2000000);
  optree->Branch("subleadci6cpfindex",&t_subleadci6cpfindex,"subleadci6cpfindex/F",2000000);
  optree->Branch("leadci6cpfmva",&t_leadci6cpfmva,"leadci6cpfmva/F",2000000);
  optree->Branch("leadci6cpfmvacat",&t_leadci6cpfmvacat,"leadci6cpfmvacat/F",2000000);
  optree->Branch("leadci6cpfmvaptom",&t_leadci6cpfmvaptom,"leadci6cpfmvaptom/F",2000000);
  optree->Branch("leadci6cpfmvaptom2",&t_leadci6cpfmvaptom2,"leadci6cpfmvaptom2/F",2000000);
  optree->Branch("subleadci6cpfmva",&t_subleadci6cpfmva,"subleadci6cpfmva/F",2000000);
  optree->Branch("subleadci6cpfmvacat",&t_subleadci6cpfmvacat,"subleadci6cpfmvacat/F",2000000);
  optree->Branch("subleadci6cpfmvaptom",&t_subleadci6cpfmvaptom,"subleadci6cpfmvaptom/F",2000000);
  optree->Branch("subleadci6cpfmvaptom2",&t_subleadci6cpfmvaptom2,"subleadci6cpfmvaptom2/F",2000000);
  optree->Branch("diphomva",&t_diphomva,"diphomva/F",2000000);
  optree->Branch("diphomva2",&t_diphomva2,"diphomva2/F",2000000);
  optree->Branch("leadci4cindex",&t_leadci4cindex,"leadci4cindex/F",2000000);
  optree->Branch("subleadci4cindex",&t_subleadci4cindex,"subleadci4cindex/F",2000000);

  /*     jgb
  optree->Branch("lead_beta_recvtx_out30",&t_lead_beta_recvtx_out30,"lead_beta_recvtx_out30/F",2000000);
  optree->Branch("lead_beta_recvtx_out40",&t_lead_beta_recvtx_out40,"lead_beta_recvtx_out40/F",2000000);
  optree->Branch("lead_beta_badvtx_out30",&t_lead_beta_badvtx_out30,"lead_beta_badvtx_out30/F",2000000);
  optree->Branch("lead_beta_badvtx_out40",&t_lead_beta_badvtx_out40,"lead_beta_badvtx_out40/F",2000000);
  optree->Branch("sublead_beta_recvtx_out30",&t_sublead_beta_recvtx_out30,"sublead_beta_recvtx_out30/F",2000000);
  optree->Branch("sublead_beta_recvtx_out40",&t_sublead_beta_recvtx_out40,"sublead_beta_recvtx_out40/F",2000000);
  optree->Branch("sublead_beta_badvtx_out30",&t_sublead_beta_badvtx_out30,"sublead_beta_badvtx_out30/F",2000000);
  optree->Branch("sublead_beta_badvtx_out40",&t_sublead_beta_badvtx_out40,"sublead_beta_badvtx_out40/F",2000000);
  */

  /*     jgb
  optree->Branch("lead_ecal_in035out20",&t_lead_ecal_in035out20,"lead_ecal_in035out20/F",2000000);
  optree->Branch("lead_ecal_in035out25",&t_lead_ecal_in035out25,"lead_ecal_in035out25/F",2000000);
  optree->Branch("lead_ecal_in035out30",&t_lead_ecal_in035out30,"lead_ecal_in035out30/F",2000000);
  optree->Branch("lead_ecal_in035out35",&t_lead_ecal_in035out35,"lead_ecal_in035out35/F",2000000);
  optree->Branch("lead_ecal_in035out40",&t_lead_ecal_in035out40,"lead_ecal_in035out40/F",2000000);
  optree->Branch("lead_ecal_in035out45",&t_lead_ecal_in035out45,"lead_ecal_in035out45/F",2000000);
  optree->Branch("lead_ecal_in035out50",&t_lead_ecal_in035out50,"lead_ecal_in035out50/F",2000000);
  optree->Branch("sublead_ecal_in035out20",&t_sublead_ecal_in035out20,"sublead_ecal_in035out20/F",2000000);
  optree->Branch("sublead_ecal_in035out25",&t_sublead_ecal_in035out25,"sublead_ecal_in035out25/F",2000000);
  optree->Branch("sublead_ecal_in035out30",&t_sublead_ecal_in035out30,"sublead_ecal_in035out30/F",2000000);
  optree->Branch("sublead_ecal_in035out35",&t_sublead_ecal_in035out35,"sublead_ecal_in035out35/F",2000000);
  optree->Branch("sublead_ecal_in035out40",&t_sublead_ecal_in035out40,"sublead_ecal_in035out40/F",2000000);
  optree->Branch("sublead_ecal_in035out45",&t_sublead_ecal_in035out45,"sublead_ecal_in035out45/F",2000000);
  optree->Branch("sublead_ecal_in035out50",&t_sublead_ecal_in035out50,"sublead_ecal_in035out50/F",2000000);
  */

  /*        jgb
  optree->Branch("lead_ecal_rhocorr12_in035out20",&t_lead_ecal_rhocorr12_in035out20,"lead_ecal_rhocorr12_in035out20/F",2000000);
  optree->Branch("lead_ecal_rhocorr12_in035out25",&t_lead_ecal_rhocorr12_in035out25,"lead_ecal_rhocorr12_in035out25/F",2000000);
  optree->Branch("lead_ecal_rhocorr12_in035out30",&t_lead_ecal_rhocorr12_in035out30,"lead_ecal_rhocorr12_in035out30/F",2000000);
  optree->Branch("lead_ecal_rhocorr12_in035out35",&t_lead_ecal_rhocorr12_in035out35,"lead_ecal_rhocorr12_in035out35/F",2000000);
  optree->Branch("lead_ecal_rhocorr12_in035out40",&t_lead_ecal_rhocorr12_in035out40,"lead_ecal_rhocorr12_in035out40/F",2000000);
  optree->Branch("lead_ecal_rhocorr12_in035out45",&t_lead_ecal_rhocorr12_in035out45,"lead_ecal_rhocorr12_in035out45/F",2000000);
  optree->Branch("lead_ecal_rhocorr12_in035out50",&t_lead_ecal_rhocorr12_in035out50,"lead_ecal_rhocorr12_in035out50/F",2000000);
  optree->Branch("sublead_ecal_rhocorr12_in035out20",&t_sublead_ecal_rhocorr12_in035out20,"sublead_ecal_rhocorr12_in035out20/F",2000000);
  optree->Branch("sublead_ecal_rhocorr12_in035out25",&t_sublead_ecal_rhocorr12_in035out25,"sublead_ecal_rhocorr12_in035out25/F",2000000);
  optree->Branch("sublead_ecal_rhocorr12_in035out30",&t_sublead_ecal_rhocorr12_in035out30,"sublead_ecal_rhocorr12_in035out30/F",2000000);
  optree->Branch("sublead_ecal_rhocorr12_in035out35",&t_sublead_ecal_rhocorr12_in035out35,"sublead_ecal_rhocorr12_in035out35/F",2000000);
  optree->Branch("sublead_ecal_rhocorr12_in035out40",&t_sublead_ecal_rhocorr12_in035out40,"sublead_ecal_rhocorr12_in035out40/F",2000000);
  optree->Branch("sublead_ecal_rhocorr12_in035out45",&t_sublead_ecal_rhocorr12_in035out45,"sublead_ecal_rhocorr12_in035out45/F",2000000);
  optree->Branch("sublead_ecal_rhocorr12_in035out50",&t_sublead_ecal_rhocorr12_in035out50,"sublead_ecal_rhocorr12_in035out50/F",2000000);
  */

  //sublead other
  optree->Branch("subleadtrkplusecal",&t_subleadtrkplusecal,"subleadtrkplusecal/F",2000000);

  optree->Branch("leadpixel",&t_leadpixel,"leadpixel/I",2000000);
  optree->Branch("leadsieie",&t_leadsieie,"leadsieie/F",2000000);
  optree->Branch("leadtrkhollowdr03",&t_leadtrkhollowdr03,"leadtrkhollowdr03/F",2000000);
  optree->Branch("leadtrkhollowdr04",&t_leadtrkhollowdr04,"leadtrkhollowdr04/F",2000000);
  optree->Branch("leadtrksoliddr03",&t_leadtrksoliddr03,"leadtrksoliddr03/F",2000000);
  optree->Branch("leadtrksoliddr04",&t_leadtrksoliddr04,"leadtrksoliddr04/F",2000000);
  optree->Branch("leadecaldr03",&t_leadecaldr03,"leadecaldr03/F",2000000);
  optree->Branch("leadecaldr04",&t_leadecaldr04,"leadecaldr04/F",2000000);
  optree->Branch("leadhcaldr03",&t_leadhcaldr03,"leadhcaldr03/F",2000000);
  optree->Branch("leadhcaldr04",&t_leadhcaldr04,"leadhcaldr04/F",2000000);
  optree->Branch("subleadpixel",&t_subleadpixel,"subleadpixel/I",2000000);
  optree->Branch("subleadsieie",&t_subleadsieie,"subleadsieie/F",2000000);
  optree->Branch("subleadtrkhollowdr03",&t_subleadtrkhollowdr03,"subleadtrkhollowdr03/F",2000000);
  optree->Branch("subleadtrkhollowdr04",&t_subleadtrkhollowdr04,"subleadtrkhollowdr04/F",2000000);
  optree->Branch("subleadtrksoliddr03",&t_subleadtrksoliddr03,"subleadtrksoliddr03/F",2000000);
  optree->Branch("subleadtrksoliddr04",&t_subleadtrksoliddr04,"subleadtrksoliddr04/F",2000000);
  optree->Branch("subleadecaldr03",&t_subleadecaldr03,"subleadecaldr03/F",2000000);
  optree->Branch("subleadecaldr04",&t_subleadecaldr04,"subleadecaldr04/F",2000000);
  optree->Branch("subleadhcaldr03",&t_subleadhcaldr03,"subleadhcaldr03/F",2000000);
  optree->Branch("subleadhcaldr04",&t_subleadhcaldr04,"subleadhcaldr04/F",2000000);


  // default reproduced
  /*           jgb
  optree->Branch("lead_tkiso_recvtx_030_004_0015_02_01",&t_lead_tkiso_recvtx_030_004_0015_02_01,"lead_tkiso_recvtx_030_004_0015_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0015_02_01",&t_sublead_tkiso_recvtx_030_004_0015_02_01,"sublead_tkiso_recvtx_030_004_0015_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_000_0015_02_01",&t_lead_tkiso_recvtx_030_000_0015_02_01,"lead_tkiso_recvtx_030_000_0015_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_000_0015_02_01",&t_sublead_tkiso_recvtx_030_000_0015_02_01,"sublead_tkiso_recvtx_030_000_0015_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_004_0015_02_01",&t_lead_tkiso_recvtx_040_004_0015_02_01,"lead_tkiso_recvtx_040_004_0015_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_004_0015_02_01",&t_sublead_tkiso_recvtx_040_004_0015_02_01,"sublead_tkiso_recvtx_040_004_0015_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_000_0015_02_01",&t_lead_tkiso_recvtx_040_000_0015_02_01,"lead_tkiso_recvtx_040_000_0015_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_000_0015_02_01",&t_sublead_tkiso_recvtx_040_000_0015_02_01,"sublead_tkiso_recvtx_040_000_0015_02_01/F",2000000);
  */

  //outer cone
  /*                          jgb
  optree->Branch("lead_tkiso_recvtx_026_004_0000_02_01",&t_lead_tkiso_recvtx_026_004_0000_02_01,"lead_tkiso_recvtx_026_004_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_026_004_0000_02_01",&t_sublead_tkiso_recvtx_026_004_0000_02_01,"sublead_tkiso_recvtx_026_004_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_028_004_0000_02_01",&t_lead_tkiso_recvtx_028_004_0000_02_01,"lead_tkiso_recvtx_028_004_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_028_004_0000_02_01",&t_sublead_tkiso_recvtx_028_004_0000_02_01,"sublead_tkiso_recvtx_028_004_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_032_004_0000_02_01",&t_lead_tkiso_recvtx_032_004_0000_02_01,"lead_tkiso_recvtx_032_004_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_032_004_0000_02_01",&t_sublead_tkiso_recvtx_032_004_0000_02_01,"sublead_tkiso_recvtx_032_004_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_034_004_0000_02_01",&t_lead_tkiso_recvtx_034_004_0000_02_01,"lead_tkiso_recvtx_034_004_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_034_004_0000_02_01",&t_sublead_tkiso_recvtx_034_004_0000_02_01,"sublead_tkiso_recvtx_034_004_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_036_004_0000_02_01",&t_lead_tkiso_recvtx_036_004_0000_02_01,"lead_tkiso_recvtx_036_004_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_036_004_0000_02_01",&t_sublead_tkiso_recvtx_036_004_0000_02_01,"sublead_tkiso_recvtx_036_004_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_026_004_0000_10_01",&t_lead_tkiso_recvtx_026_004_0000_10_01,"lead_tkiso_recvtx_026_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_026_004_0000_10_01",&t_sublead_tkiso_recvtx_026_004_0000_10_01,"sublead_tkiso_recvtx_026_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_028_004_0000_10_01",&t_lead_tkiso_recvtx_028_004_0000_10_01,"lead_tkiso_recvtx_028_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_028_004_0000_10_01",&t_sublead_tkiso_recvtx_028_004_0000_10_01,"sublead_tkiso_recvtx_028_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_032_004_0000_10_01",&t_lead_tkiso_recvtx_032_004_0000_10_01,"lead_tkiso_recvtx_032_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_032_004_0000_10_01",&t_sublead_tkiso_recvtx_032_004_0000_10_01,"sublead_tkiso_recvtx_032_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_034_004_0000_10_01",&t_lead_tkiso_recvtx_034_004_0000_10_01,"lead_tkiso_recvtx_034_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_034_004_0000_10_01",&t_sublead_tkiso_recvtx_034_004_0000_10_01,"sublead_tkiso_recvtx_034_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_036_004_0000_10_01",&t_lead_tkiso_recvtx_036_004_0000_10_01,"lead_tkiso_recvtx_036_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_036_004_0000_10_01",&t_sublead_tkiso_recvtx_036_004_0000_10_01,"sublead_tkiso_recvtx_036_004_0000_10_01/F",2000000);
  */

  //inner cone
  /*         jgb
  optree->Branch("lead_tkiso_recvtx_030_003_0000_02_01",&t_lead_tkiso_recvtx_030_003_0000_02_01,"lead_tkiso_recvtx_030_003_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_003_0000_02_01",&t_sublead_tkiso_recvtx_030_003_0000_02_01,"sublead_tkiso_recvtx_030_003_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_005_0000_02_01",&t_lead_tkiso_recvtx_030_005_0000_02_01,"lead_tkiso_recvtx_030_005_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_005_0000_02_01",&t_sublead_tkiso_recvtx_030_005_0000_02_01,"sublead_tkiso_recvtx_030_005_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_006_0000_02_01",&t_lead_tkiso_recvtx_030_006_0000_02_01,"lead_tkiso_recvtx_030_006_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_006_0000_02_01",&t_sublead_tkiso_recvtx_030_006_0000_02_01,"sublead_tkiso_recvtx_030_006_0000_02_01/F",2000000);
  */

  /*        jgb
  optree->Branch("lead_tkiso_recvtx_030_000_0000_10_01",&t_lead_tkiso_recvtx_030_000_0000_10_01,"lead_tkiso_recvtx_030_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_001_0000_10_01",&t_lead_tkiso_recvtx_030_001_0000_10_01,"lead_tkiso_recvtx_030_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_002_0000_10_01",&t_lead_tkiso_recvtx_030_002_0000_10_01,"lead_tkiso_recvtx_030_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_003_0000_10_01",&t_lead_tkiso_recvtx_030_003_0000_10_01,"lead_tkiso_recvtx_030_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_035_000_0000_10_01",&t_lead_tkiso_recvtx_035_000_0000_10_01,"lead_tkiso_recvtx_035_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_035_001_0000_10_01",&t_lead_tkiso_recvtx_035_001_0000_10_01,"lead_tkiso_recvtx_035_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_035_002_0000_10_01",&t_lead_tkiso_recvtx_035_002_0000_10_01,"lead_tkiso_recvtx_035_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_035_003_0000_10_01",&t_lead_tkiso_recvtx_035_003_0000_10_01,"lead_tkiso_recvtx_035_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_000_0000_10_01",&t_lead_tkiso_recvtx_040_000_0000_10_01,"lead_tkiso_recvtx_040_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_001_0000_10_01",&t_lead_tkiso_recvtx_040_001_0000_10_01,"lead_tkiso_recvtx_040_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_002_0000_10_01",&t_lead_tkiso_recvtx_040_002_0000_10_01,"lead_tkiso_recvtx_040_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_003_0000_10_01",&t_lead_tkiso_recvtx_040_003_0000_10_01,"lead_tkiso_recvtx_040_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_000_0000_10_01",&t_sublead_tkiso_recvtx_030_000_0000_10_01,"sublead_tkiso_recvtx_030_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_001_0000_10_01",&t_sublead_tkiso_recvtx_030_001_0000_10_01,"sublead_tkiso_recvtx_030_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_002_0000_10_01",&t_sublead_tkiso_recvtx_030_002_0000_10_01,"sublead_tkiso_recvtx_030_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_003_0000_10_01",&t_sublead_tkiso_recvtx_030_003_0000_10_01,"sublead_tkiso_recvtx_030_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_035_000_0000_10_01",&t_sublead_tkiso_recvtx_035_000_0000_10_01,"sublead_tkiso_recvtx_035_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_035_001_0000_10_01",&t_sublead_tkiso_recvtx_035_001_0000_10_01,"sublead_tkiso_recvtx_035_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_035_002_0000_10_01",&t_sublead_tkiso_recvtx_035_002_0000_10_01,"sublead_tkiso_recvtx_035_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_035_003_0000_10_01",&t_sublead_tkiso_recvtx_035_003_0000_10_01,"sublead_tkiso_recvtx_035_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_000_0000_10_01",&t_sublead_tkiso_recvtx_040_000_0000_10_01,"sublead_tkiso_recvtx_040_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_001_0000_10_01",&t_sublead_tkiso_recvtx_040_001_0000_10_01,"sublead_tkiso_recvtx_040_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_002_0000_10_01",&t_sublead_tkiso_recvtx_040_002_0000_10_01,"sublead_tkiso_recvtx_040_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_003_0000_10_01",&t_sublead_tkiso_recvtx_040_003_0000_10_01,"sublead_tkiso_recvtx_040_003_0000_10_01/F",2000000);
  */

  /*       jgb
  optree->Branch("lead_tkiso_genvtx_030_000_0000_10_01",&t_lead_tkiso_genvtx_030_000_0000_10_01,"lead_tkiso_genvtx_030_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_030_001_0000_10_01",&t_lead_tkiso_genvtx_030_001_0000_10_01,"lead_tkiso_genvtx_030_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_030_002_0000_10_01",&t_lead_tkiso_genvtx_030_002_0000_10_01,"lead_tkiso_genvtx_030_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_030_003_0000_10_01",&t_lead_tkiso_genvtx_030_003_0000_10_01,"lead_tkiso_genvtx_030_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_035_000_0000_10_01",&t_lead_tkiso_genvtx_035_000_0000_10_01,"lead_tkiso_genvtx_035_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_035_001_0000_10_01",&t_lead_tkiso_genvtx_035_001_0000_10_01,"lead_tkiso_genvtx_035_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_035_002_0000_10_01",&t_lead_tkiso_genvtx_035_002_0000_10_01,"lead_tkiso_genvtx_035_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_035_003_0000_10_01",&t_lead_tkiso_genvtx_035_003_0000_10_01,"lead_tkiso_genvtx_035_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_040_000_0000_10_01",&t_lead_tkiso_genvtx_040_000_0000_10_01,"lead_tkiso_genvtx_040_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_040_001_0000_10_01",&t_lead_tkiso_genvtx_040_001_0000_10_01,"lead_tkiso_genvtx_040_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_040_002_0000_10_01",&t_lead_tkiso_genvtx_040_002_0000_10_01,"lead_tkiso_genvtx_040_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_genvtx_040_003_0000_10_01",&t_lead_tkiso_genvtx_040_003_0000_10_01,"lead_tkiso_genvtx_040_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_030_000_0000_10_01",&t_sublead_tkiso_genvtx_030_000_0000_10_01,"sublead_tkiso_genvtx_030_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_030_001_0000_10_01",&t_sublead_tkiso_genvtx_030_001_0000_10_01,"sublead_tkiso_genvtx_030_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_030_002_0000_10_01",&t_sublead_tkiso_genvtx_030_002_0000_10_01,"sublead_tkiso_genvtx_030_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_030_003_0000_10_01",&t_sublead_tkiso_genvtx_030_003_0000_10_01,"sublead_tkiso_genvtx_030_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_035_000_0000_10_01",&t_sublead_tkiso_genvtx_035_000_0000_10_01,"sublead_tkiso_genvtx_035_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_035_001_0000_10_01",&t_sublead_tkiso_genvtx_035_001_0000_10_01,"sublead_tkiso_genvtx_035_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_035_002_0000_10_01",&t_sublead_tkiso_genvtx_035_002_0000_10_01,"sublead_tkiso_genvtx_035_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_035_003_0000_10_01",&t_sublead_tkiso_genvtx_035_003_0000_10_01,"sublead_tkiso_genvtx_035_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_040_000_0000_10_01",&t_sublead_tkiso_genvtx_040_000_0000_10_01,"sublead_tkiso_genvtx_040_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_040_001_0000_10_01",&t_sublead_tkiso_genvtx_040_001_0000_10_01,"sublead_tkiso_genvtx_040_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_040_002_0000_10_01",&t_sublead_tkiso_genvtx_040_002_0000_10_01,"sublead_tkiso_genvtx_040_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_genvtx_040_003_0000_10_01",&t_sublead_tkiso_genvtx_040_003_0000_10_01,"sublead_tkiso_genvtx_040_003_0000_10_01/F",2000000);
  */

  /*         jgb
  optree->Branch("lead_tkiso_badvtx_030_000_0000_10_01",&t_lead_tkiso_badvtx_030_000_0000_10_01,"lead_tkiso_badvtx_030_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_001_0000_10_01",&t_lead_tkiso_badvtx_030_001_0000_10_01,"lead_tkiso_badvtx_030_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_002_0000_10_01",&t_lead_tkiso_badvtx_030_002_0000_10_01,"lead_tkiso_badvtx_030_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_003_0000_10_01",&t_lead_tkiso_badvtx_030_003_0000_10_01,"lead_tkiso_badvtx_030_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_035_000_0000_10_01",&t_lead_tkiso_badvtx_035_000_0000_10_01,"lead_tkiso_badvtx_035_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_035_001_0000_10_01",&t_lead_tkiso_badvtx_035_001_0000_10_01,"lead_tkiso_badvtx_035_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_035_002_0000_10_01",&t_lead_tkiso_badvtx_035_002_0000_10_01,"lead_tkiso_badvtx_035_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_035_003_0000_10_01",&t_lead_tkiso_badvtx_035_003_0000_10_01,"lead_tkiso_badvtx_035_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_000_0000_10_01",&t_lead_tkiso_badvtx_040_000_0000_10_01,"lead_tkiso_badvtx_040_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_001_0000_10_01",&t_lead_tkiso_badvtx_040_001_0000_10_01,"lead_tkiso_badvtx_040_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_002_0000_10_01",&t_lead_tkiso_badvtx_040_002_0000_10_01,"lead_tkiso_badvtx_040_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_003_0000_10_01",&t_lead_tkiso_badvtx_040_003_0000_10_01,"lead_tkiso_badvtx_040_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_000_0000_10_01",&t_sublead_tkiso_badvtx_030_000_0000_10_01,"sublead_tkiso_badvtx_030_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_001_0000_10_01",&t_sublead_tkiso_badvtx_030_001_0000_10_01,"sublead_tkiso_badvtx_030_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_002_0000_10_01",&t_sublead_tkiso_badvtx_030_002_0000_10_01,"sublead_tkiso_badvtx_030_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_003_0000_10_01",&t_sublead_tkiso_badvtx_030_003_0000_10_01,"sublead_tkiso_badvtx_030_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_035_000_0000_10_01",&t_sublead_tkiso_badvtx_035_000_0000_10_01,"sublead_tkiso_badvtx_035_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_035_001_0000_10_01",&t_sublead_tkiso_badvtx_035_001_0000_10_01,"sublead_tkiso_badvtx_035_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_035_002_0000_10_01",&t_sublead_tkiso_badvtx_035_002_0000_10_01,"sublead_tkiso_badvtx_035_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_035_003_0000_10_01",&t_sublead_tkiso_badvtx_035_003_0000_10_01,"sublead_tkiso_badvtx_035_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_000_0000_10_01",&t_sublead_tkiso_badvtx_040_000_0000_10_01,"sublead_tkiso_badvtx_040_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_001_0000_10_01",&t_sublead_tkiso_badvtx_040_001_0000_10_01,"sublead_tkiso_badvtx_040_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_002_0000_10_01",&t_sublead_tkiso_badvtx_040_002_0000_10_01,"sublead_tkiso_badvtx_040_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_003_0000_10_01",&t_sublead_tkiso_badvtx_040_003_0000_10_01,"sublead_tkiso_badvtx_040_003_0000_10_01/F",2000000);
  */

  /*         jgb
  optree->Branch("lead_tkiso_badvtx_rhocorr314_030_000_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_030_000_0000_10_01,"lead_tkiso_badvtx_rhocorr314_030_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_030_001_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_030_001_0000_10_01,"lead_tkiso_badvtx_rhocorr314_030_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_030_002_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_030_002_0000_10_01,"lead_tkiso_badvtx_rhocorr314_030_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_030_003_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_030_003_0000_10_01,"lead_tkiso_badvtx_rhocorr314_030_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_035_000_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_035_000_0000_10_01,"lead_tkiso_badvtx_rhocorr314_035_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_035_001_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_035_001_0000_10_01,"lead_tkiso_badvtx_rhocorr314_035_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_035_002_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_035_002_0000_10_01,"lead_tkiso_badvtx_rhocorr314_035_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_035_003_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_035_003_0000_10_01,"lead_tkiso_badvtx_rhocorr314_035_003_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_040_000_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_040_000_0000_10_01,"lead_tkiso_badvtx_rhocorr314_040_000_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_040_001_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_040_001_0000_10_01,"lead_tkiso_badvtx_rhocorr314_040_001_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_040_002_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_040_002_0000_10_01,"lead_tkiso_badvtx_rhocorr314_040_002_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_rhocorr314_040_003_0000_10_01",&t_lead_tkiso_badvtx_rhocorr314_040_003_0000_10_01,"lead_tkiso_badvtx_rhocorr314_040_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_030_000_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_030_000_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_030_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_030_001_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_030_001_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_030_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_030_002_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_030_002_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_030_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_030_003_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_030_003_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_030_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_035_000_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_035_000_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_035_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_035_001_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_035_001_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_035_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_035_002_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_035_002_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_035_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_035_003_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_035_003_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_035_003_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_040_000_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_040_000_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_040_000_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_040_001_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_040_001_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_040_001_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_040_002_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_040_002_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_040_002_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_rhocorr314_040_003_0000_10_01",&t_sublead_tkiso_badvtx_rhocorr314_040_003_0000_10_01,"sublead_tkiso_badvtx_rhocorr314_040_003_0000_10_01/F",2000000);
  */

  optree->Branch("lead_tkiso_recvtx_030_005_0000_10_01",&t_lead_tkiso_recvtx_030_005_0000_10_01,"lead_tkiso_recvtx_030_005_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_005_0000_10_01",&t_sublead_tkiso_recvtx_030_005_0000_10_01,"sublead_tkiso_recvtx_030_005_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_006_0000_10_01",&t_lead_tkiso_recvtx_030_006_0000_10_01,"lead_tkiso_recvtx_030_006_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_006_0000_10_01",&t_sublead_tkiso_recvtx_030_006_0000_10_01,"sublead_tkiso_recvtx_030_006_0000_10_01/F",2000000);

  //eta strip
  /*        jgb
  optree->Branch("lead_tkiso_recvtx_030_004_0010_02_01",&t_lead_tkiso_recvtx_030_004_0010_02_01,"lead_tkiso_recvtx_030_004_0010_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0010_02_01",&t_sublead_tkiso_recvtx_030_004_0010_02_01,"sublead_tkiso_recvtx_030_004_0010_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_004_0020_02_01",&t_lead_tkiso_recvtx_030_004_0020_02_01,"lead_tkiso_recvtx_030_004_0020_02_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0020_02_01",&t_sublead_tkiso_recvtx_030_004_0020_02_01,"sublead_tkiso_recvtx_030_004_0020_02_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_004_0010_10_01",&t_lead_tkiso_recvtx_030_004_0010_10_01,"lead_tkiso_recvtx_030_004_0010_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0010_10_01",&t_sublead_tkiso_recvtx_030_004_0010_10_01,"sublead_tkiso_recvtx_030_004_0010_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_004_0020_10_01",&t_lead_tkiso_recvtx_030_004_0020_10_01,"lead_tkiso_recvtx_030_004_0020_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0020_10_01",&t_sublead_tkiso_recvtx_030_004_0020_10_01,"sublead_tkiso_recvtx_030_004_0020_10_01/F",2000000);

  // deltaz
  optree->Branch("lead_tkiso_recvtx_030_004_0000_06_01",&t_lead_tkiso_recvtx_030_004_0000_06_01,"lead_tkiso_recvtx_030_004_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0000_06_01",&t_sublead_tkiso_recvtx_030_004_0000_06_01,"sublead_tkiso_recvtx_030_004_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_004_0000_10_01",&t_lead_tkiso_recvtx_030_004_0000_10_01,"lead_tkiso_recvtx_030_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0000_10_01",&t_sublead_tkiso_recvtx_030_004_0000_10_01,"sublead_tkiso_recvtx_030_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_004_0000_15_01",&t_lead_tkiso_recvtx_030_004_0000_15_01,"lead_tkiso_recvtx_030_004_0000_15_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0000_15_01",&t_sublead_tkiso_recvtx_030_004_0000_15_01,"sublead_tkiso_recvtx_030_004_0000_15_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_004_0000_99_01",&t_lead_tkiso_recvtx_030_004_0000_99_01,"lead_tkiso_recvtx_030_004_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0000_99_01",&t_sublead_tkiso_recvtx_030_004_0000_99_01,"sublead_tkiso_recvtx_030_004_0000_99_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_004_0000_06_01",&t_lead_tkiso_recvtx_040_004_0000_06_01,"lead_tkiso_recvtx_040_004_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_004_0000_06_01",&t_sublead_tkiso_recvtx_040_004_0000_06_01,"sublead_tkiso_recvtx_040_004_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_004_0000_10_01",&t_lead_tkiso_recvtx_040_004_0000_10_01,"lead_tkiso_recvtx_040_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_004_0000_10_01",&t_sublead_tkiso_recvtx_040_004_0000_10_01,"sublead_tkiso_recvtx_040_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_004_0000_99_01",&t_lead_tkiso_recvtx_040_004_0000_99_01,"lead_tkiso_recvtx_040_004_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_004_0000_99_01",&t_sublead_tkiso_recvtx_040_004_0000_99_01,"sublead_tkiso_recvtx_040_004_0000_99_01/F",2000000);

  optree->Branch("lead_tkiso_badvtx_030_004_0000_02_01",&t_lead_tkiso_badvtx_030_004_0000_02_01,"lead_tkiso_badvtx_030_004_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_004_0000_02_01",&t_sublead_tkiso_badvtx_030_004_0000_02_01,"sublead_tkiso_badvtx_030_004_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_004_0000_06_01",&t_lead_tkiso_badvtx_030_004_0000_06_01,"lead_tkiso_badvtx_030_004_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_004_0000_06_01",&t_sublead_tkiso_badvtx_030_004_0000_06_01,"sublead_tkiso_badvtx_030_004_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_004_0000_10_01",&t_lead_tkiso_badvtx_030_004_0000_10_01,"lead_tkiso_badvtx_030_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_004_0000_10_01",&t_sublead_tkiso_badvtx_030_004_0000_10_01,"sublead_tkiso_badvtx_030_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_004_0000_15_01",&t_lead_tkiso_badvtx_030_004_0000_15_01,"lead_tkiso_badvtx_030_004_0000_15_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_004_0000_15_01",&t_sublead_tkiso_badvtx_030_004_0000_15_01,"sublead_tkiso_badvtx_030_004_0000_15_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_004_0000_99_01",&t_lead_tkiso_badvtx_030_004_0000_99_01,"lead_tkiso_badvtx_030_004_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_004_0000_99_01",&t_sublead_tkiso_badvtx_030_004_0000_99_01,"sublead_tkiso_badvtx_030_004_0000_99_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_004_0000_02_01",&t_lead_tkiso_badvtx_040_004_0000_02_01,"lead_tkiso_badvtx_040_004_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_004_0000_02_01",&t_sublead_tkiso_badvtx_040_004_0000_02_01,"sublead_tkiso_badvtx_040_004_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_004_0000_06_01",&t_lead_tkiso_badvtx_040_004_0000_06_01,"lead_tkiso_badvtx_040_004_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_004_0000_06_01",&t_sublead_tkiso_badvtx_040_004_0000_06_01,"sublead_tkiso_badvtx_040_004_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_004_0000_10_01",&t_lead_tkiso_badvtx_040_004_0000_10_01,"lead_tkiso_badvtx_040_004_0000_10_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_004_0000_10_01",&t_sublead_tkiso_badvtx_040_004_0000_10_01,"sublead_tkiso_badvtx_040_004_0000_10_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_004_0000_99_01",&t_lead_tkiso_badvtx_040_004_0000_99_01,"lead_tkiso_badvtx_040_004_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_004_0000_99_01",&t_sublead_tkiso_badvtx_040_004_0000_99_01,"sublead_tkiso_badvtx_040_004_0000_99_01/F",2000000);

  optree->Branch("lead_tkiso_recvtx_030_000_0000_06_01",&t_lead_tkiso_recvtx_030_000_0000_06_01,"lead_tkiso_recvtx_030_000_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_000_0000_06_01",&t_sublead_tkiso_recvtx_030_000_0000_06_01,"sublead_tkiso_recvtx_030_000_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_030_000_0000_99_01",&t_lead_tkiso_recvtx_030_000_0000_99_01,"lead_tkiso_recvtx_030_000_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_000_0000_99_01",&t_sublead_tkiso_recvtx_030_000_0000_99_01,"sublead_tkiso_recvtx_030_000_0000_99_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_000_0000_06_01",&t_lead_tkiso_recvtx_040_000_0000_06_01,"lead_tkiso_recvtx_040_000_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_000_0000_06_01",&t_sublead_tkiso_recvtx_040_000_0000_06_01,"sublead_tkiso_recvtx_040_000_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_recvtx_040_000_0000_99_01",&t_lead_tkiso_recvtx_040_000_0000_99_01,"lead_tkiso_recvtx_040_000_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_040_000_0000_99_01",&t_sublead_tkiso_recvtx_040_000_0000_99_01,"sublead_tkiso_recvtx_040_000_0000_99_01/F",2000000);

  optree->Branch("lead_tkiso_badvtx_030_000_0000_02_01",&t_lead_tkiso_badvtx_030_000_0000_02_01,"lead_tkiso_badvtx_030_000_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_000_0000_02_01",&t_sublead_tkiso_badvtx_030_000_0000_02_01,"sublead_tkiso_badvtx_030_000_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_000_0000_06_01",&t_lead_tkiso_badvtx_030_000_0000_06_01,"lead_tkiso_badvtx_030_000_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_000_0000_06_01",&t_sublead_tkiso_badvtx_030_000_0000_06_01,"sublead_tkiso_badvtx_030_000_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_030_000_0000_99_01",&t_lead_tkiso_badvtx_030_000_0000_99_01,"lead_tkiso_badvtx_030_000_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_030_000_0000_99_01",&t_sublead_tkiso_badvtx_030_000_0000_99_01,"sublead_tkiso_badvtx_030_000_0000_99_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_000_0000_02_01",&t_lead_tkiso_badvtx_040_000_0000_02_01,"lead_tkiso_badvtx_040_000_0000_02_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_000_0000_02_01",&t_sublead_tkiso_badvtx_040_000_0000_02_01,"sublead_tkiso_badvtx_040_000_0000_02_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_000_0000_06_01",&t_lead_tkiso_badvtx_040_000_0000_06_01,"lead_tkiso_badvtx_040_000_0000_06_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_000_0000_06_01",&t_sublead_tkiso_badvtx_040_000_0000_06_01,"sublead_tkiso_badvtx_040_000_0000_06_01/F",2000000);
  optree->Branch("lead_tkiso_badvtx_040_000_0000_99_01",&t_lead_tkiso_badvtx_040_000_0000_99_01,"lead_tkiso_badvtx_040_000_0000_99_01/F",2000000);
  optree->Branch("sublead_tkiso_badvtx_040_000_0000_99_01",&t_sublead_tkiso_badvtx_040_000_0000_99_01,"sublead_tkiso_badvtx_040_000_0000_99_01/F",2000000);

  // dxy
  optree->Branch("lead_tkiso_recvtx_030_004_0000_10_02",&t_lead_tkiso_recvtx_030_004_0000_10_02,"lead_tkiso_recvtx_030_004_0000_10_02/F",2000000);
  optree->Branch("sublead_tkiso_recvtx_030_004_0000_10_02",&t_sublead_tkiso_recvtx_030_004_0000_10_02,"sublead_tkiso_recvtx_030_004_0000_10_02/F",2000000);
  */

  // vertex z, sumpt
  optree->Branch("genvtxz",&t_genvtxz,"genvtxz/F",2000000);
  optree->Branch("recvtxz",&t_recvtxz,"recvtxz/F",2000000);
  optree->Branch("badvtxz",&t_badvtxz,"badvtxz/F",2000000);
  optree->Branch("genvtx_sumpt",&t_genvtx_sumpt,"genvtx_sumpt/F",2000000);
  optree->Branch("recvtx_sumpt",&t_recvtx_sumpt,"recvtx_sumpt/F",2000000);
  optree->Branch("badvtx_sumpt",&t_badvtx_sumpt,"badvtx_sumpt/F",2000000);

  // pfiso
  optree->Branch("lead_pfiso_charged03",&t_lead_pfiso_charged03,"lead_pfiso_charged03/F",2000000);
  optree->Branch("lead_pfiso_photon03",&t_lead_pfiso_photon03,"lead_pfiso_photon03/F",2000000);
  optree->Branch("lead_pfiso_neutral03",&t_lead_pfiso_neutral03,"lead_pfiso_neutral03/F",2000000);
  optree->Branch("lead_pfiso_charged04",&t_lead_pfiso_charged04,"lead_pfiso_charged04/F",2000000);
  optree->Branch("lead_pfiso_photon04",&t_lead_pfiso_photon04,"lead_pfiso_photon04/F",2000000);
  optree->Branch("lead_pfiso_neutral04",&t_lead_pfiso_neutral04,"lead_pfiso_neutral04/F",2000000);
  optree->Branch("sublead_pfiso_charged03",&t_sublead_pfiso_charged03,"sublead_pfiso_charged03/F",2000000);
  optree->Branch("sublead_pfiso_photon03",&t_sublead_pfiso_photon03,"sublead_pfiso_photon03/F",2000000);
  optree->Branch("sublead_pfiso_neutral03",&t_sublead_pfiso_neutral03,"sublead_pfiso_neutral03/F",2000000);
  optree->Branch("sublead_pfiso_charged04",&t_sublead_pfiso_charged04,"sublead_pfiso_charged04/F",2000000);
  optree->Branch("sublead_pfiso_photon04",&t_sublead_pfiso_photon04,"sublead_pfiso_photon04/F",2000000);
  optree->Branch("sublead_pfiso_neutral04",&t_sublead_pfiso_neutral04,"sublead_pfiso_neutral04/F",2000000);

  optree->Branch("lead_pfiso_charged_badvtx_03",&t_lead_pfiso_charged_badvtx_03,"lead_pfiso_charged_badvtx_03/F",2000000);
  optree->Branch("lead_pfiso_photon_badvtx_03",&t_lead_pfiso_photon_badvtx_03,"lead_pfiso_photon_badvtx_03/F",2000000);
  optree->Branch("lead_pfiso_neutral_badvtx_03",&t_lead_pfiso_neutral_badvtx_03,"lead_pfiso_neutral_badvtx_03/F",2000000);
  optree->Branch("lead_pfiso_charged_badvtx_04",&t_lead_pfiso_charged_badvtx_04,"lead_pfiso_charged_badvtx_04/F",2000000);
  optree->Branch("lead_pfiso_photon_badvtx_04",&t_lead_pfiso_photon_badvtx_04,"lead_pfiso_photon_badvtx_04/F",2000000);
  optree->Branch("lead_pfiso_neutral_badvtx_04",&t_lead_pfiso_neutral_badvtx_04,"lead_pfiso_neutral_badvtx_04/F",2000000);
  optree->Branch("sublead_pfiso_charged_badvtx_03",&t_sublead_pfiso_charged_badvtx_03,"sublead_pfiso_charged_badvtx_03/F",2000000);
  optree->Branch("sublead_pfiso_photon_badvtx_03",&t_sublead_pfiso_photon_badvtx_03,"sublead_pfiso_photon_badvtx_03/F",2000000);
  optree->Branch("sublead_pfiso_neutral_badvtx_03",&t_sublead_pfiso_neutral_badvtx_03,"sublead_pfiso_neutral_badvtx_03/F",2000000);
  optree->Branch("sublead_pfiso_charged_badvtx_04",&t_sublead_pfiso_charged_badvtx_04,"sublead_pfiso_charged_badvtx_04/F",2000000);
  optree->Branch("sublead_pfiso_photon_badvtx_04",&t_sublead_pfiso_photon_badvtx_04,"sublead_pfiso_photon_badvtx_04/F",2000000);
  optree->Branch("sublead_pfiso_neutral_badvtx_04",&t_sublead_pfiso_neutral_badvtx_04,"sublead_pfiso_neutral_badvtx_04/F",2000000);


  //delta-r to track
  /*         jgb
  optree->Branch("lead_drtotk_10_06",&t_lead_drtotk_10_06,"lead_drtotk_10_06/F",2000000);
  optree->Branch("sublead_drtotk_10_06",&t_sublead_drtotk_10_06,"sublead_drtotk_10_06/F",2000000);
  optree->Branch("lead_drtotk_15_06",&t_lead_drtotk_15_06,"lead_drtotk_15_06/F",2000000);
  optree->Branch("sublead_drtotk_15_06",&t_sublead_drtotk_15_06,"sublead_drtotk_15_06/F",2000000);
  optree->Branch("lead_drtotk_20_06",&t_lead_drtotk_20_06,"lead_drtotk_20_06/F",2000000);
  optree->Branch("sublead_drtotk_20_06",&t_sublead_drtotk_20_06,"sublead_drtotk_20_06/F",2000000);
  optree->Branch("lead_drtotk_25_06",&t_lead_drtotk_25_06,"lead_drtotk_25_06/F",2000000);
  optree->Branch("sublead_drtotk_25_06",&t_sublead_drtotk_25_06,"sublead_drtotk_25_06/F",2000000);
  optree->Branch("lead_drtotkworstvtx_10_06",&t_lead_drtotkworstvtx_10_06,"lead_drtotkworstvtx_10_06/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_10_06",&t_sublead_drtotkworstvtx_10_06,"sublead_drtotkworstvtx_10_06/F",2000000);
  optree->Branch("lead_drtotkworstvtx_15_06",&t_lead_drtotkworstvtx_15_06,"lead_drtotkworstvtx_15_06/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_15_06",&t_sublead_drtotkworstvtx_15_06,"sublead_drtotkworstvtx_15_06/F",2000000);
  optree->Branch("lead_drtotkworstvtx_20_06",&t_lead_drtotkworstvtx_20_06,"lead_drtotkworstvtx_20_06/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_20_06",&t_sublead_drtotkworstvtx_20_06,"sublead_drtotkworstvtx_20_06/F",2000000);
  optree->Branch("lead_drtotkworstvtx_25_06",&t_lead_drtotkworstvtx_25_06,"lead_drtotkworstvtx_25_06/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_25_06",&t_sublead_drtotkworstvtx_25_06,"sublead_drtotkworstvtx_25_06/F",2000000);

  optree->Branch("lead_drtotk_10_10",&t_lead_drtotk_10_10,"lead_drtotk_10_10/F",2000000);
  optree->Branch("sublead_drtotk_10_10",&t_sublead_drtotk_10_10,"sublead_drtotk_10_10/F",2000000);
  optree->Branch("lead_drtotk_15_10",&t_lead_drtotk_15_10,"lead_drtotk_15_10/F",2000000);
  optree->Branch("sublead_drtotk_15_10",&t_sublead_drtotk_15_10,"sublead_drtotk_15_10/F",2000000);
  optree->Branch("lead_drtotk_20_10",&t_lead_drtotk_20_10,"lead_drtotk_20_10/F",2000000);
  optree->Branch("sublead_drtotk_20_10",&t_sublead_drtotk_20_10,"sublead_drtotk_20_10/F",2000000);
  optree->Branch("lead_drtotk_25_10",&t_lead_drtotk_25_10,"lead_drtotk_25_10/F",2000000);
  optree->Branch("sublead_drtotk_25_10",&t_sublead_drtotk_25_10,"sublead_drtotk_25_10/F",2000000);
  optree->Branch("lead_drtotkworstvtx_10_10",&t_lead_drtotkworstvtx_10_10,"lead_drtotkworstvtx_10_10/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_10_10",&t_sublead_drtotkworstvtx_10_10,"sublead_drtotkworstvtx_10_10/F",2000000);
  optree->Branch("lead_drtotkworstvtx_15_10",&t_lead_drtotkworstvtx_15_10,"lead_drtotkworstvtx_15_10/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_15_10",&t_sublead_drtotkworstvtx_15_10,"sublead_drtotkworstvtx_15_10/F",2000000);
  optree->Branch("lead_drtotkworstvtx_20_10",&t_lead_drtotkworstvtx_20_10,"lead_drtotkworstvtx_20_10/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_20_10",&t_sublead_drtotkworstvtx_20_10,"sublead_drtotkworstvtx_20_10/F",2000000);
  optree->Branch("lead_drtotkworstvtx_25_10",&t_lead_drtotkworstvtx_25_10,"lead_drtotkworstvtx_25_10/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_25_10",&t_sublead_drtotkworstvtx_25_10,"sublead_drtotkworstvtx_25_10/F",2000000);

  optree->Branch("lead_drtotk_10_99",&t_lead_drtotk_10_99,"lead_drtotk_10_99/F",2000000);
  optree->Branch("sublead_drtotk_10_99",&t_sublead_drtotk_10_99,"sublead_drtotk_10_99/F",2000000);
  optree->Branch("lead_drtotk_15_99",&t_lead_drtotk_15_99,"lead_drtotk_15_99/F",2000000);
  optree->Branch("sublead_drtotk_15_99",&t_sublead_drtotk_15_99,"sublead_drtotk_15_99/F",2000000);
  optree->Branch("lead_drtotk_20_99",&t_lead_drtotk_20_99,"lead_drtotk_20_99/F",2000000);
  optree->Branch("sublead_drtotk_20_99",&t_sublead_drtotk_20_99,"sublead_drtotk_20_99/F",2000000);
  optree->Branch("lead_drtotkworstvtx_10_99",&t_lead_drtotkworstvtx_10_99,"lead_drtotkworstvtx_10_99/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_10_99",&t_sublead_drtotkworstvtx_10_99,"sublead_drtotkworstvtx_10_99/F",2000000);
  optree->Branch("lead_drtotkworstvtx_15_99",&t_lead_drtotkworstvtx_15_99,"lead_drtotkworstvtx_15_99/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_15_99",&t_sublead_drtotkworstvtx_15_99,"sublead_drtotkworstvtx_15_99/F",2000000);
  optree->Branch("lead_drtotkworstvtx_20_99",&t_lead_drtotkworstvtx_20_99,"lead_drtotkworstvtx_20_99/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_20_99",&t_sublead_drtotkworstvtx_20_99,"sublead_drtotkworstvtx_20_99/F",2000000);
  optree->Branch("lead_drtotkworstvtx_25_99",&t_lead_drtotkworstvtx_25_99,"lead_drtotkworstvtx_25_99/F",2000000);
  optree->Branch("sublead_drtotkworstvtx_25_99",&t_sublead_drtotkworstvtx_25_99,"sublead_drtotkworstvtx_25_99/F",2000000);
  */
  optree->Branch("lead_drtotk_25_99",&t_lead_drtotk_25_99,"lead_drtotk_25_99/F",2000000);
  optree->Branch("sublead_drtotk_25_99",&t_sublead_drtotk_25_99,"sublead_drtotk_25_99/F",2000000);

  /*     
    optree->Branch("lpt_n", &lpt_n, "lpt_n/I",2000000);
    optree->Branch("lpt_emu_n", &lpt_emu_n, "lpt_emu_n/I",2000000);
    optree->Branch("lpt_mu_n", &lpt_mu_n, "lpt_mu_n/I",2000000);
    optree->Branch("lpt_el_n", &lpt_el_n, "lpt_el_n/I",2000000);
    optree->Branch("lpt_pho_n", &lpt_pho_n, "lpt_pho_n/I",2000000);
    optree->Branch("lpt_pdgid", &lpt_pdgid, "lpt_pdgid[lpt_n]/I",2000000);
    optree->Branch("lpt_ind", &lpt_ind, "lpt_ind[lpt_n]/I",2000000);
    optree->Branch("lpt_duplicate", &lpt_duplicate, "lpt_duplicate[lpt_n]/I",2000000);
    optree->Branch("lpt_indgen", &lpt_indgen, "lpt_indgen[lpt_n]/I",2000000);
    optree->Branch("lpt_drmatch", &lpt_drmatch, "lpt_drmatch[lpt_n]/F",2000000);
    lpt_p4  = new TClonesArray("TLorentzVector", MAX_LEPTONS);
    optree->Branch("lpt_p4", "TClonesArray", &lpt_p4, 10000000, 2);

    optree->Branch("met_pfmet", &met_pfmet, "met_pfmet/F",2000000); 
    optree->Branch("met_phi_pfmet", &met_phi_pfmet, "met_phi_pfmet/F",2000000); 
  */

}



// Jim's Functions BEGIN

// Jim's Functions BEGIN

void StatAnalysisExclusive::setDiphoCuts() {
  //  Sets the diphoton cut levels (from CiC diphoton cuts) just called once.
  //  The output of the cut-setting program can just be copied here to update the numbers. 
  //  the new cuts are after the k-factors greatly reducing the high PT higgs

 //----------------------------sob value=0.000444444     end of iteration 9          UltraLoose  ---------------------------
             float cutdiphoptom0[] = {  0.0022,   0.00147,     0.004,    0.0034};
                  float cutdmom0[] = {   0.033,     0.029,     0.026,     0.027};
      float cutfsubleadcutindex0[] = {     8.5,       8.5,       8.5,       8.5};
         float cutfleadcutindex0[] = {     8.5,       7.5,       8.5,       6.5};
                float cutetamax0[] = {    1.44,      1.44,       2.5,       2.5};
        float cutsubleadptomass0[] = {    0.25,      0.25,      0.25,      0.25};
               float cutsumptom0[] = {     0.6,      0.61,      0.59,       0.6};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.945325        fake=0.729655        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.000777778     end of iteration 7          VeryLoose  ---------------------------
             float cutdiphoptom1[] = {  0.0022,   0.00173,    0.0052,    0.0058};
                  float cutdmom1[] = {   0.033,     0.022,     0.026,     0.024};
      float cutfsubleadcutindex1[] = {     8.5,       8.5,       8.5,       8.5};
         float cutfleadcutindex1[] = {     7.5,       6.5,       6.5,       5.5};
                float cutetamax1[] = {    1.44,      1.44,       2.5,       2.4};
        float cutsubleadptomass1[] = {    0.25,      0.25,      0.25,      0.25};
               float cutsumptom1[] = {     0.6,      0.61,       0.6,      0.62};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.868998        fake=0.547487        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00122222     end of iteration 7          Loose  ---------------------------
             float cutdiphoptom2[] = {  0.0022,   0.00173,    0.0052,    0.0191};
                  float cutdmom2[] = {   0.033,    0.0194,     0.026,     0.022};
      float cutfsubleadcutindex2[] = {     7.5,       7.5,       7.5,       7.5};
         float cutfleadcutindex2[] = {     5.5,       5.5,       4.5,       4.5};
                float cutetamax2[] = {    1.44,      1.44,       2.4,       2.4};
        float cutsubleadptomass2[] = {    0.25,      0.25,      0.25,      0.25};
               float cutsumptom2[] = {     0.6,      0.65,      0.61,      0.74};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.728761        fake=0.347994        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00155556     end of iteration 7          Medium  ---------------------------
             float cutdiphoptom3[] = {  0.0022,   0.00197,    0.0052,     0.026};
                  float cutdmom3[] = {  0.0164,    0.0174,     0.026,     0.022};
      float cutfsubleadcutindex3[] = {     7.5,       7.5,       7.5,       7.5};
         float cutfleadcutindex3[] = {     5.5,       5.5,       3.5,       4.5};
                float cutetamax3[] = {    1.44,      1.44,       2.3,       2.2};
        float cutsubleadptomass3[] = {    0.25,      0.25,      0.26,      0.31};
               float cutsumptom3[] = {    0.61,       0.7,      0.72,      0.74};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.684399        fake=0.298078        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00211111     end of iteration 7          Tight  ---------------------------
             float cutdiphoptom4[] = {  0.0025,    0.0023,    0.0052,     0.145};
                  float cutdmom4[] = {  0.0141,    0.0162,     0.021,     0.021};
      float cutfsubleadcutindex4[] = {     5.5,       6.5,       6.5,       7.5};
         float cutfleadcutindex4[] = {     5.5,       4.5,       3.5,       3.5};
                float cutetamax4[] = {    1.44,      1.44,       2.3,       2.2};
        float cutsubleadptomass4[] = {    0.25,      0.25,      0.33,      0.35};
               float cutsumptom4[] = {    0.62,       0.7,      0.73,       0.8};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.553555        fake=0.188396        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00266667     end of iteration 7          SuperTight  ---------------------------
             float cutdiphoptom5[] = {  0.0025,    0.0042,    0.0143,      0.22};
                  float cutdmom5[] = {   0.014,    0.0156,      0.02,    0.0198};
      float cutfsubleadcutindex5[] = {     5.5,       6.5,       5.5,       7.5};
         float cutfleadcutindex5[] = {     3.5,       3.5,       3.5,       3.5};
                float cutetamax5[] = {    1.44,      1.44,       2.2,       2.2};
        float cutsubleadptomass5[] = {    0.25,      0.25,      0.39,      0.38};
               float cutsumptom5[] = {    0.72,      0.73,      0.79,       0.8};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.442634        fake=0.122665        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00333333     end of iteration 7          HyperTight0  ---------------------------
             float cutdiphoptom6[] = {  0.0028,     0.021,     0.029,      0.35};
                  float cutdmom6[] = {  0.0137,    0.0146,      0.02,    0.0198};
      float cutfsubleadcutindex6[] = {     5.5,       6.5,       3.5,       7.5};
         float cutfleadcutindex6[] = {     2.5,       3.5,       3.5,       3.5};
                float cutetamax6[] = {    1.44,      1.44,       2.2,       2.2};
        float cutsubleadptomass6[] = {    0.25,      0.25,      0.43,      0.42};
               float cutsumptom6[] = {    0.72,      0.83,      0.89,      0.97};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.370691        fake=0.0871436        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00388889     end of iteration 7          HyperTight1  ---------------------------
             float cutdiphoptom7[] = {  0.0028,     0.053,     0.073,      0.35};
                  float cutdmom7[] = {  0.0128,    0.0141,      0.02,    0.0195};
      float cutfsubleadcutindex7[] = {     5.5,       6.5,       3.5,       7.5};
         float cutfleadcutindex7[] = {     1.5,       3.5,       3.5,       3.5};
                float cutetamax7[] = {    1.44,      1.44,       2.2,       2.1};
        float cutsubleadptomass7[] = {    0.25,      0.25,      0.43,      0.42};
               float cutsumptom7[] = {     0.8,      0.84,      0.89,      0.99};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.336263        fake=0.072843        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00466667     end of iteration 7          HyperTight2  ---------------------------
             float cutdiphoptom8[] = {  0.0049,     0.067,     0.073,      0.35};
                  float cutdmom8[] = {  0.0118,    0.0141,      0.02,     0.019};
      float cutfsubleadcutindex8[] = {     5.5,       4.5,      -0.5,       7.5};
         float cutfleadcutindex8[] = {     1.5,       3.5,       3.5,       3.5};
                float cutetamax8[] = {    1.44,      1.44,       2.2,       2.1};
        float cutsubleadptomass8[] = {    0.25,      0.26,      0.43,      0.42};
               float cutsumptom8[] = {    0.85,      0.88,      0.89,         1};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.269197        fake=0.0507364        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00544444     end of iteration 7          HyperTight3  ---------------------------
             float cutdiphoptom9[] = {  0.0049,     0.081,     0.073,      0.35};
                  float cutdmom9[] = {  0.0118,    0.0132,      0.02,    0.0189};
      float cutfsubleadcutindex9[] = {     4.5,       4.5,      -0.5,       7.5};
         float cutfleadcutindex9[] = {     1.5,       3.5,       3.5,       2.5};
                float cutetamax9[] = {    1.44,      1.44,       2.2,       2.1};
        float cutsubleadptomass9[] = {    0.25,      0.26,      0.43,      0.43};
               float cutsumptom9[] = {    0.87,      0.88,      0.89,         1};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.22419        fake=0.0373714        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00633333     end of iteration 7          HyperTight4  ---------------------------
            float cutdiphoptom10[] = {   0.008,     0.085,     0.073,      0.35};
                 float cutdmom10[] = {  0.0116,    0.0132,      0.02,    0.0187};
     float cutfsubleadcutindex10[] = {     4.5,       4.5,      -0.5,       7.5};
        float cutfleadcutindex10[] = {     1.5,       1.5,      -0.5,       2.5};
               float cutetamax10[] = {    1.44,      1.44,       2.2,       2.1};
       float cutsubleadptomass10[] = {    0.25,      0.27,      0.43,      0.47};
              float cutsumptom10[] = {    0.89,      0.88,      0.89,         1};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.183736        fake=0.0283085        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.00733333     end of iteration 7          HyperTight5  ---------------------------
            float cutdiphoptom11[] = {   0.044,     0.128,     0.073,      0.35};
                 float cutdmom11[] = {  0.0116,    0.0132,      0.02,    0.0167};
     float cutfsubleadcutindex11[] = {     4.5,       4.5,      -0.5,       4.5};
        float cutfleadcutindex11[] = {     1.5,       1.5,      -0.5,      -0.5};
               float cutetamax11[] = {    1.44,      1.44,       2.2,      1.94};
       float cutsubleadptomass11[] = {    0.29,      0.38,      0.43,       0.5};
              float cutsumptom11[] = {    0.92,      0.95,      0.89,         1};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.106966        fake=0.012933        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.008     end of iteration 7          HyperTight6  ---------------------------
            float cutdiphoptom12[] = {   0.051,     0.192,     0.073,      0.35};
                 float cutdmom12[] = {  0.0116,    0.0132,      0.02,    0.0167};
     float cutfsubleadcutindex12[] = {     4.5,       4.5,      -0.5,       2.5};
        float cutfleadcutindex12[] = {     1.5,       1.5,      -0.5,      -0.5};
               float cutetamax12[] = {    1.44,      1.44,       2.2,      1.94};
       float cutsubleadptomass12[] = {     0.3,       0.4,      0.43,       0.5};
              float cutsumptom12[] = {    0.94,      0.99,      0.89,         1};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0919332        fake=0.0100011        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //----------------------------sob value=0.0122222     end of iteration 7          HyperTight7  ---------------------------
            float cutdiphoptom13[] = {     0.2,      0.37,     0.073,      0.35};
                 float cutdmom13[] = {  0.0113,    0.0132,      0.02,    0.0167};
     float cutfsubleadcutindex13[] = {     2.5,       4.5,      -0.5,       2.5};
        float cutfleadcutindex13[] = {     1.5,       1.5,      -0.5,      -0.5};
               float cutetamax13[] = {    1.44,      1.44,       2.2,      1.94};
       float cutsubleadptomass13[] = {    0.41,      0.41,      0.43,       0.5};
              float cutsumptom13[] = {    0.94,      1.07,      0.89,         1};
 //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0427014        fake=0.00253804        <<<<<<<<<<<<<<<<<<<<<<<<<<<<

  for(int ic=0;ic!=4;++ic) {
    cutdiphoptom[0][ic]=cutdiphoptom0[ic];    cutdiphoptom[1][ic]=cutdiphoptom1[ic];    cutdiphoptom[2][ic]=cutdiphoptom2[ic];    cutdiphoptom[3][ic]=cutdiphoptom3[ic];    cutdiphoptom[4][ic]=cutdiphoptom4[ic];  
    cutdiphoptom[5][ic]=cutdiphoptom5[ic];    cutdiphoptom[6][ic]=cutdiphoptom6[ic];    cutdiphoptom[7][ic]=cutdiphoptom7[ic];    cutdiphoptom[8][ic]=cutdiphoptom8[ic];    cutdiphoptom[9][ic]=cutdiphoptom9[ic];  
    cutdiphoptom[10][ic]=cutdiphoptom10[ic];  cutdiphoptom[11][ic]=cutdiphoptom11[ic];  cutdiphoptom[12][ic]=cutdiphoptom12[ic];  cutdiphoptom[13][ic]=cutdiphoptom13[ic]; 

    cutdmom[0][ic]=cutdmom0[ic];    cutdmom[1][ic]=cutdmom1[ic];    cutdmom[2][ic]=cutdmom2[ic];    cutdmom[3][ic]=cutdmom3[ic];    cutdmom[4][ic]=cutdmom4[ic];  
    cutdmom[5][ic]=cutdmom5[ic];    cutdmom[6][ic]=cutdmom6[ic];    cutdmom[7][ic]=cutdmom7[ic];    cutdmom[8][ic]=cutdmom8[ic];    cutdmom[9][ic]=cutdmom9[ic];  
    cutdmom[10][ic]=cutdmom10[ic];  cutdmom[11][ic]=cutdmom11[ic];  cutdmom[12][ic]=cutdmom12[ic];  cutdmom[13][ic]=cutdmom13[ic]; 

    cutfsubleadcutindex[0][ic]=cutfsubleadcutindex0[ic];    cutfsubleadcutindex[1][ic]=cutfsubleadcutindex1[ic];    cutfsubleadcutindex[2][ic]=cutfsubleadcutindex2[ic];    cutfsubleadcutindex[3][ic]=cutfsubleadcutindex3[ic];    cutfsubleadcutindex[4][ic]=cutfsubleadcutindex4[ic];  
    cutfsubleadcutindex[5][ic]=cutfsubleadcutindex5[ic];    cutfsubleadcutindex[6][ic]=cutfsubleadcutindex6[ic];    cutfsubleadcutindex[7][ic]=cutfsubleadcutindex7[ic];    cutfsubleadcutindex[8][ic]=cutfsubleadcutindex8[ic];    cutfsubleadcutindex[9][ic]=cutfsubleadcutindex9[ic];  
    cutfsubleadcutindex[10][ic]=cutfsubleadcutindex10[ic];  cutfsubleadcutindex[11][ic]=cutfsubleadcutindex11[ic];  cutfsubleadcutindex[12][ic]=cutfsubleadcutindex12[ic];   cutfsubleadcutindex[13][ic]=cutfsubleadcutindex13[ic]; 

    cutfleadcutindex[0][ic]=cutfleadcutindex0[ic];    cutfleadcutindex[1][ic]=cutfleadcutindex1[ic];    cutfleadcutindex[2][ic]=cutfleadcutindex2[ic];    cutfleadcutindex[3][ic]=cutfleadcutindex3[ic];    cutfleadcutindex[4][ic]=cutfleadcutindex4[ic];  
    cutfleadcutindex[5][ic]=cutfleadcutindex5[ic];    cutfleadcutindex[6][ic]=cutfleadcutindex6[ic];    cutfleadcutindex[7][ic]=cutfleadcutindex7[ic];    cutfleadcutindex[8][ic]=cutfleadcutindex8[ic];    cutfleadcutindex[9][ic]=cutfleadcutindex9[ic];  
    cutfleadcutindex[10][ic]=cutfleadcutindex10[ic];  cutfleadcutindex[11][ic]=cutfleadcutindex11[ic];  cutfleadcutindex[12][ic]=cutfleadcutindex12[ic];   cutfleadcutindex[13][ic]=cutfleadcutindex13[ic]; 

    cutetamax[0][ic]=cutetamax0[ic];    cutetamax[1][ic]=cutetamax1[ic];    cutetamax[2][ic]=cutetamax2[ic];    cutetamax[3][ic]=cutetamax3[ic];    cutetamax[4][ic]=cutetamax4[ic];  
    cutetamax[5][ic]=cutetamax5[ic];    cutetamax[6][ic]=cutetamax6[ic];    cutetamax[7][ic]=cutetamax7[ic];    cutetamax[8][ic]=cutetamax8[ic];    cutetamax[9][ic]=cutetamax9[ic];  
    cutetamax[10][ic]=cutetamax10[ic];  cutetamax[11][ic]=cutetamax11[ic];  cutetamax[12][ic]=cutetamax12[ic];  cutetamax[13][ic]=cutetamax13[ic]; 

    cutsubleadptomass[0][ic]=cutsubleadptomass0[ic];    cutsubleadptomass[1][ic]=cutsubleadptomass1[ic];    cutsubleadptomass[2][ic]=cutsubleadptomass2[ic];    cutsubleadptomass[3][ic]=cutsubleadptomass3[ic];    cutsubleadptomass[4][ic]=cutsubleadptomass4[ic];  
    cutsubleadptomass[5][ic]=cutsubleadptomass5[ic];    cutsubleadptomass[6][ic]=cutsubleadptomass6[ic];    cutsubleadptomass[7][ic]=cutsubleadptomass7[ic];    cutsubleadptomass[8][ic]=cutsubleadptomass8[ic];    cutsubleadptomass[9][ic]=cutsubleadptomass9[ic];  
    cutsubleadptomass[10][ic]=cutsubleadptomass10[ic];  cutsubleadptomass[11][ic]=cutsubleadptomass11[ic];  cutsubleadptomass[12][ic]=cutsubleadptomass12[ic];  cutsubleadptomass[13][ic]=cutsubleadptomass13[ic];

    cutsumptom[0][ic]=cutsumptom0[ic];    cutsumptom[1][ic]=cutsumptom1[ic];    cutsumptom[2][ic]=cutsumptom2[ic];    cutsumptom[3][ic]=cutsumptom3[ic];    cutsumptom[4][ic]=cutsumptom4[ic];  
    cutsumptom[5][ic]=cutsumptom5[ic];    cutsumptom[6][ic]=cutsumptom6[ic];    cutsumptom[7][ic]=cutsumptom7[ic];    cutsumptom[8][ic]=cutsumptom8[ic];    cutsumptom[9][ic]=cutsumptom9[ic];  
    cutsumptom[10][ic]=cutsumptom10[ic];  cutsumptom[11][ic]=cutsumptom11[ic];  cutsumptom[12][ic]=cutsumptom12[ic];  cutsumptom[13][ic]=cutsumptom13[ic]; 

  }

}


void StatAnalysisExclusive::setPhotonCuts6() {
  //  Inits the cut levels in 6 categories for photon selection.  Just called once.
  //  Output of cut-setting program can just be copied here.

//  after marco's changes to the diphoton selection there are less events on output tree:  retune
//  ----------------------------sob value=0.000254777     end of iteration 9          Loose  ---------------------------
      float cutsubleadisosumoet0[] = {      22,        23,       9.3,      19.8,         6,       4.9};
   float cutsubleadisosumoetbad0[] = {      78,        87,      13.8,        85,      14.6,      10.7};
    float cutsubleadtrkisooetom0[] = {    14.1,        14,       5.9,      12.8,       3.2,       2.7};
          float cutsubleadsieie0[] = {  0.0117,    0.0129,    0.0122,     0.032,     0.033,     0.034};
         float cutsubleadhovere0[] = {   0.114,      0.14,     0.137,     0.146,     0.114,     0.126};
             float cutsubleadr90[] = {    0.94,       0.9,      0.23,      0.94,       0.9,      0.23};
  float cutsublead_drtotk_25_990[] = {    0.99,      0.99,    0.0105,    0.0086,     0.009,    0.0092};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.967005        fake=0.717552        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00044586     end of iteration 8          Medium  ---------------------------
      float cutsubleadisosumoet1[] = {     9.8,      16.7,       6.6,       6.5,       5.2,       4.5};
   float cutsubleadisosumoetbad1[] = {      59,        87,      11.8,        24,       9.5,       7.5};
    float cutsubleadtrkisooetom1[] = {     6.7,       9.4,       3.9,       4.6,       3.1,      1.96};
          float cutsubleadsieie1[] = {  0.0114,    0.0126,    0.0122,     0.031,     0.033,     0.029};
         float cutsubleadhovere1[] = {   0.107,      0.14,     0.137,     0.146,     0.108,      0.06};
             float cutsubleadr91[] = {    0.94,       0.9,      0.23,      0.94,       0.9,      0.23};
  float cutsublead_drtotk_25_991[] = {       1,         1,    0.0128,    0.0168,    0.0138,    0.0182};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.938526        fake=0.620558        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.000636943     end of iteration 8          Tight  ---------------------------
      float cutsubleadisosumoet2[] = {     7.9,       7.3,       5.7,       5.7,       4.9,       4.5};
   float cutsubleadisosumoetbad2[] = {      54,      16.5,       9.8,        22,       8.5,       6.8};
    float cutsubleadtrkisooetom2[] = {     4.9,       6.4,       2.9,         4,       3.1,      1.91};
          float cutsubleadsieie2[] = {  0.0114,    0.0126,    0.0112,     0.028,     0.033,     0.029};
         float cutsubleadhovere2[] = {   0.106,      0.14,     0.105,     0.146,     0.108,      0.06};
             float cutsubleadr92[] = {    0.94,       0.9,      0.23,      0.94,       0.9,      0.23};
  float cutsublead_drtotk_25_992[] = {       1,         1,    0.0153,      0.99,      0.67,      0.32};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.901715        fake=0.53066        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.000828025     end of iteration 8          SuperTight  ---------------------------
      float cutsubleadisosumoet3[] = {     7.2,       6.9,       5.6,       4.6,       4.9,       4.5};
   float cutsubleadisosumoetbad3[] = {      20,      11.5,       8.8,        13,       7.5,       6.1};
    float cutsubleadtrkisooetom3[] = {     4.6,       3.8,       2.9,       3.1,       2.3,      1.91};
          float cutsubleadsieie3[] = {   0.011,    0.0124,    0.0108,     0.028,     0.028,     0.029};
         float cutsubleadhovere3[] = {   0.104,      0.14,     0.105,     0.146,     0.081,     0.059};
             float cutsubleadr93[] = {    0.94,       0.9,      0.24,      0.94,       0.9,      0.27};
  float cutsublead_drtotk_25_993[] = {       1,         1,     0.154,         1,         1,      0.82};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.867332        fake=0.465838        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.0010828     end of iteration 8          HyperTight1  ---------------------------
      float cutsubleadisosumoet4[] = {     6.5,       6.9,       5.5,       3.9,         4,       4.4};
   float cutsubleadisosumoetbad4[] = {      20,        10,       7.6,       9.1,       6.1,       5.4};
    float cutsubleadtrkisooetom4[] = {     3.2,       3.8,       2.9,       2.6,       2.3,      1.86};
          float cutsubleadsieie4[] = {  0.0109,    0.0108,    0.0106,     0.028,     0.028,     0.029};
         float cutsubleadhovere4[] = {   0.096,     0.139,     0.105,     0.146,     0.081,     0.041};
             float cutsubleadr94[] = {    0.94,       0.9,      0.24,      0.94,       0.9,      0.39};
  float cutsublead_drtotk_25_994[] = {       1,         1,     0.154,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.822623        fake=0.407654        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00140127     end of iteration 8          HyperTight2  ---------------------------
      float cutsubleadisosumoet5[] = {     6.1,       6.9,       5.5,       3.7,       3.7,       3.3};
   float cutsubleadisosumoetbad5[] = {    10.9,       9.3,       6.9,       5.8,       5.4,       4.8};
    float cutsubleadtrkisooetom5[] = {       3,       3.8,       2.9,       2.1,      1.42,      0.98};
          float cutsubleadsieie5[] = {  0.0108,    0.0105,    0.0101,     0.028,     0.027,     0.028};
         float cutsubleadhovere5[] = {   0.086,     0.139,      0.07,     0.146,     0.081,     0.039};
             float cutsubleadr95[] = {    0.94,       0.9,      0.25,      0.94,       0.9,      0.39};
  float cutsublead_drtotk_25_995[] = {       1,         1,      0.21,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.754229        fake=0.337185        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00171975     end of iteration 8          HyperTight3  ---------------------------
      float cutsubleadisosumoet6[] = {     5.9,       5.8,         5,       3.4,       2.6,       3.3};
   float cutsubleadisosumoetbad6[] = {     9.6,       8.7,       6.5,       5.2,       4.8,       4.3};
    float cutsubleadtrkisooetom6[] = {     2.5,       3.6,       2.5,       2.1,      0.82,      0.86};
          float cutsubleadsieie6[] = {  0.0108,    0.0105,    0.0101,     0.028,     0.027,     0.028};
         float cutsubleadhovere6[] = {   0.064,     0.139,      0.07,     0.116,     0.041,     0.022};
             float cutsubleadr96[] = {    0.94,       0.9,      0.25,      0.94,       0.9,      0.39};
  float cutsublead_drtotk_25_996[] = {       1,         1,      0.99,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.690855        fake=0.282639        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.0022293     end of iteration 8          HyperTight4  ---------------------------
      float cutsubleadisosumoet7[] = {       5,       5.3,       4.6,       2.8,       2.1,       2.8};
   float cutsubleadisosumoetbad7[] = {     8.1,       7.1,       6.1,       5.2,       3.8,         4};
    float cutsubleadtrkisooetom7[] = {     2.5,       2.8,       2.5,      0.83,      0.58,      0.76};
          float cutsubleadsieie7[] = {  0.0106,    0.0105,    0.0099,     0.028,     0.026,     0.028};
         float cutsubleadhovere7[] = {   0.045,     0.063,     0.038,     0.021,     0.028,    0.0183};
             float cutsubleadr97[] = {    0.94,       0.9,      0.26,      0.94,       0.9,      0.52};
  float cutsublead_drtotk_25_997[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.557679        fake=0.197394        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00286624     end of iteration 8          HyperTight5  ---------------------------
      float cutsubleadisosumoet8[] = {     4.7,       4.6,       3.8,       2.4,      1.88,       2.6};
   float cutsubleadisosumoetbad8[] = {     6.4,       6.2,       5.6,       3.3,       3.3,       3.5};
    float cutsubleadtrkisooetom8[] = {     2.1,       2.5,      1.71,      0.31,      0.47,       0.7};
          float cutsubleadsieie8[] = {  0.0106,    0.0104,    0.0099,     0.027,     0.026,     0.027};
         float cutsubleadhovere8[] = {   0.045,     0.062,      0.03,   0.00021,     0.027,    0.0105};
             float cutsubleadr98[] = {    0.94,       0.9,      0.26,      0.95,       0.9,      0.81};
  float cutsublead_drtotk_25_998[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.427211        fake=0.128498        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00382166     end of iteration 8          HyperTight6  ---------------------------
      float cutsubleadisosumoet9[] = {     4.3,       3.8,       3.8,      1.76,      1.54,       2.4};
   float cutsubleadisosumoetbad9[] = {     5.6,       5.8,       4.9,         3,       3.3,       3.5};
    float cutsubleadtrkisooetom9[] = {    1.95,       2.5,      1.58,      0.31,     0.142,      0.41};
          float cutsubleadsieie9[] = {  0.0104,    0.0103,    0.0099,     0.026,     0.026,     0.025};
         float cutsubleadhovere9[] = {   0.045,     0.049,     0.027,  0.000129,     0.027,    0.0105};
             float cutsubleadr99[] = {    0.94,       0.9,      0.42,      0.95,       0.9,      0.85};
  float cutsublead_drtotk_25_999[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.309945        fake=0.0811574        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00477707     end of iteration 8          HyperTight7  ---------------------------
     float cutsubleadisosumoet10[] = {     3.3,       3.8,       3.7,      1.55,      1.51,       2.4};
  float cutsubleadisosumoetbad10[] = {     5.6,       5.2,       4.1,       2.4,       3.3,       3.5};
   float cutsubleadtrkisooetom10[] = {    0.89,      1.79,      1.46,      0.31,     0.142,      0.41};
         float cutsubleadsieie10[] = {  0.0103,    0.0103,    0.0098,     0.025,     0.026,     0.025};
        float cutsubleadhovere10[] = {   0.027,     0.042,     0.027,  0.000129,    0.0096,    0.0105};
            float cutsubleadr910[] = {    0.94,       0.9,      0.42,      0.95,      0.91,      0.85};
 float cutsublead_drtotk_25_9910[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.206397        fake=0.0487522        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00573248     end of iteration 8          HyperTight8  ---------------------------
     float cutsubleadisosumoet11[] = {     3.1,       3.2,       2.3,      1.55,      1.51,       2.4};
  float cutsubleadisosumoetbad11[] = {     3.5,       5.2,       3.3,       2.1,       3.3,       3.5};
   float cutsubleadtrkisooetom11[] = {    0.42,      0.62,      0.73,      0.31,     0.142,      0.41};
         float cutsubleadsieie11[] = {  0.0103,    0.0103,    0.0095,     0.025,     0.025,     0.025};
        float cutsubleadhovere11[] = {   0.027,     0.036,     0.016,   1.7e-05,    0.0096,    0.0105};
            float cutsubleadr911[] = {    0.94,       0.9,      0.47,      0.96,      0.91,      0.85};
 float cutsublead_drtotk_25_9911[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0619041        fake=0.0131031        <<<<<<<<<<<<<<<<<<<<<<<<<<<<



  for(int ic=0; ic<6; ic++) {

    cutsubleadisosumoet6c[0][ic]=cutsubleadisosumoet0[ic];     cutsubleadisosumoet6c[1][ic]=cutsubleadisosumoet1[ic];     cutsubleadisosumoet6c[2][ic]=cutsubleadisosumoet2[ic];     cutsubleadisosumoet6c[3][ic]=cutsubleadisosumoet3[ic];    
    cutsubleadisosumoet6c[4][ic]=cutsubleadisosumoet4[ic];     cutsubleadisosumoet6c[5][ic]=cutsubleadisosumoet5[ic];     cutsubleadisosumoet6c[6][ic]=cutsubleadisosumoet6[ic];     cutsubleadisosumoet6c[7][ic]=cutsubleadisosumoet7[ic];
    cutsubleadisosumoet6c[8][ic]=cutsubleadisosumoet8[ic];     cutsubleadisosumoet6c[9][ic]=cutsubleadisosumoet9[ic];     cutsubleadisosumoet6c[10][ic]=cutsubleadisosumoet10[ic];   

    cutsubleadisosumoetbad6c[0][ic]=cutsubleadisosumoetbad0[ic];     cutsubleadisosumoetbad6c[1][ic]=cutsubleadisosumoetbad1[ic];     cutsubleadisosumoetbad6c[2][ic]=cutsubleadisosumoetbad2[ic];     cutsubleadisosumoetbad6c[3][ic]=cutsubleadisosumoetbad3[ic];    
    cutsubleadisosumoetbad6c[4][ic]=cutsubleadisosumoetbad4[ic];     cutsubleadisosumoetbad6c[5][ic]=cutsubleadisosumoetbad5[ic];     cutsubleadisosumoetbad6c[6][ic]=cutsubleadisosumoetbad6[ic];     cutsubleadisosumoetbad6c[7][ic]=cutsubleadisosumoetbad7[ic];
    cutsubleadisosumoetbad6c[8][ic]=cutsubleadisosumoetbad8[ic];     cutsubleadisosumoetbad6c[9][ic]=cutsubleadisosumoetbad9[ic];     cutsubleadisosumoetbad6c[10][ic]=cutsubleadisosumoetbad10[ic];   
    
    cutsubleadtrkisooetom6c[0][ic]=cutsubleadtrkisooetom0[ic];     cutsubleadtrkisooetom6c[1][ic]=cutsubleadtrkisooetom1[ic];     cutsubleadtrkisooetom6c[2][ic]=cutsubleadtrkisooetom2[ic];     cutsubleadtrkisooetom6c[3][ic]=cutsubleadtrkisooetom3[ic];    
    cutsubleadtrkisooetom6c[4][ic]=cutsubleadtrkisooetom4[ic];     cutsubleadtrkisooetom6c[5][ic]=cutsubleadtrkisooetom5[ic];     cutsubleadtrkisooetom6c[6][ic]=cutsubleadtrkisooetom6[ic];     cutsubleadtrkisooetom6c[7][ic]=cutsubleadtrkisooetom7[ic];
    cutsubleadtrkisooetom6c[8][ic]=cutsubleadtrkisooetom8[ic];     cutsubleadtrkisooetom6c[9][ic]=cutsubleadtrkisooetom9[ic];     cutsubleadtrkisooetom6c[10][ic]=cutsubleadtrkisooetom10[ic];   
    
    cutsubleadsieie6c[0][ic]=cutsubleadsieie0[ic];     cutsubleadsieie6c[1][ic]=cutsubleadsieie1[ic];     cutsubleadsieie6c[2][ic]=cutsubleadsieie2[ic];     cutsubleadsieie6c[3][ic]=cutsubleadsieie3[ic];    
    cutsubleadsieie6c[4][ic]=cutsubleadsieie4[ic];     cutsubleadsieie6c[5][ic]=cutsubleadsieie5[ic];     cutsubleadsieie6c[6][ic]=cutsubleadsieie6[ic];     cutsubleadsieie6c[7][ic]=cutsubleadsieie7[ic];
    cutsubleadsieie6c[8][ic]=cutsubleadsieie8[ic];     cutsubleadsieie6c[9][ic]=cutsubleadsieie9[ic];     cutsubleadsieie6c[10][ic]=cutsubleadsieie10[ic];   
    
    cutsubleadhovere6c[0][ic]=cutsubleadhovere0[ic];     cutsubleadhovere6c[1][ic]=cutsubleadhovere1[ic];     cutsubleadhovere6c[2][ic]=cutsubleadhovere2[ic];     cutsubleadhovere6c[3][ic]=cutsubleadhovere3[ic];    
    cutsubleadhovere6c[4][ic]=cutsubleadhovere4[ic];     cutsubleadhovere6c[5][ic]=cutsubleadhovere5[ic];     cutsubleadhovere6c[6][ic]=cutsubleadhovere6[ic];     cutsubleadhovere6c[7][ic]=cutsubleadhovere7[ic];
    cutsubleadhovere6c[8][ic]=cutsubleadhovere8[ic];     cutsubleadhovere6c[9][ic]=cutsubleadhovere9[ic];     cutsubleadhovere6c[10][ic]=cutsubleadhovere10[ic];   
    
    cutsubleadr96c[0][ic]=cutsubleadr90[ic];     cutsubleadr96c[1][ic]=cutsubleadr91[ic];     cutsubleadr96c[2][ic]=cutsubleadr92[ic];     cutsubleadr96c[3][ic]=cutsubleadr93[ic];    
    cutsubleadr96c[4][ic]=cutsubleadr94[ic];     cutsubleadr96c[5][ic]=cutsubleadr95[ic];     cutsubleadr96c[6][ic]=cutsubleadr96[ic];     cutsubleadr96c[7][ic]=cutsubleadr97[ic];
    cutsubleadr96c[8][ic]=cutsubleadr98[ic];     cutsubleadr96c[9][ic]=cutsubleadr99[ic];     cutsubleadr96c[10][ic]=cutsubleadr910[ic];   

    cutsublead_drtotk_25_996c[0][ic]=cutsublead_drtotk_25_990[ic];   cutsublead_drtotk_25_996c[1][ic]=cutsublead_drtotk_25_991[ic];   cutsublead_drtotk_25_996c[2][ic]=cutsublead_drtotk_25_992[ic];   cutsublead_drtotk_25_996c[3][ic]=cutsublead_drtotk_25_993[ic];    
    cutsublead_drtotk_25_996c[4][ic]=cutsublead_drtotk_25_994[ic];   cutsublead_drtotk_25_996c[5][ic]=cutsublead_drtotk_25_995[ic];   cutsublead_drtotk_25_996c[6][ic]=cutsublead_drtotk_25_996[ic];   cutsublead_drtotk_25_996c[7][ic]=cutsublead_drtotk_25_997[ic];
    cutsublead_drtotk_25_996c[8][ic]=cutsublead_drtotk_25_998[ic];   cutsublead_drtotk_25_996c[9][ic]=cutsublead_drtotk_25_999[ic];   cutsublead_drtotk_25_996c[10][ic]=cutsublead_drtotk_25_9910[ic];   

  }
}



void StatAnalysisExclusive::setPhotonCuts6pf() {
  //  Inits the cut levels in 6 categories for photon selection pf.  Just called once.
  //  Output of cut-setting program can just be copied here.


//  ----------------------------sob value=0.000222222     end of iteration 9          Loose  ---------------------------
    float cutpfisosumoet0[] = {     9.2,       7.4,         6,       9.3,        10,         7};
 float cutpfisosumoetbad0[] = {      47,        55,      11.3,      19.3,        58,       7.7};
       float cutpfchisooet0[] = {     6.3,       4.8,       4.4,       4.4,       7.1,       4.2};
   float cutpfhcalisooetom0[] = {    11.6,        22,       8.7,       8.3,       8.9,       7.5};
          float cutpfsieie0[] = {  0.0116,    0.0114,    0.0103,     0.032,     0.031,      0.03};
         float cutpfhovere0[] = {   0.129,     0.143,     0.112,     0.146,     0.136,     0.147};
             float cutpfr90[] = {    0.94,       0.9,      0.23,      0.94,       0.9,      0.23};
  float cutpf_drtotk_25_990[] = {    0.99,      0.96,    0.0064,    0.0066,    0.0041,    0.0032};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.969349        fake=0.702721        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.000388889     end of iteration 8          Medium  ---------------------------
    float cutpfisosumoet1[] = {       8,       6.5,       5.1,       6.2,       6.7,       4.8};
 float cutpfisosumoetbad1[] = {      44,      11.4,       8.2,      10.1,       9.5,       6.5};
       float cutpfchisooet1[] = {     5.6,       3.7,       3.5,       3.7,       4.5,       3.4};
   float cutpfhcalisooetom1[] = {    11.6,        22,       7.8,       8.2,       6.5,       7.4};
          float cutpfsieie1[] = {  0.0116,    0.0114,    0.0101,     0.029,      0.03,     0.029};
         float cutpfhovere1[] = {   0.129,     0.142,      0.11,     0.146,     0.136,     0.075};
             float cutpfr91[] = {    0.94,       0.9,      0.23,      0.94,       0.9,      0.23};
  float cutpf_drtotk_25_991[] = {       1,         1,    0.0114,    0.0081,    0.0068,    0.0109};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.949368        fake=0.617595        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.000555556     end of iteration 8          Tight  ---------------------------
    float cutpfisosumoet2[] = {     7.2,       6.2,       4.7,       5.7,       6.7,       4.3};
 float cutpfisosumoetbad2[] = {    11.8,       8.2,       6.8,       8.4,       6.5,       5.6};
       float cutpfchisooet2[] = {     4.7,       3.7,       2.5,       3.4,       4.5,       2.3};
   float cutpfhcalisooetom2[] = {    11.6,        22,       7.8,       8.2,       6.4,       7.4};
          float cutpfsieie2[] = {  0.0113,    0.0111,    0.0101,     0.028,      0.03,     0.028};
         float cutpfhovere2[] = {   0.123,     0.142,      0.11,     0.146,     0.134,     0.061};
             float cutpfr92[] = {    0.94,       0.9,      0.24,      0.94,       0.9,      0.23};
  float cutpf_drtotk_25_992[] = {       1,         1,    0.0125,    0.0094,      0.02,    0.0145};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.928023        fake=0.559129        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.000722222     end of iteration 8          SuperTight  ---------------------------
    float cutpfisosumoet3[] = {     6.1,       5.5,       4.7,       5.6,       6.3,       3.8};
 float cutpfisosumoetbad3[] = {    10.9,       7.5,       6.2,       6.9,       5.6,       5.1};
       float cutpfchisooet3[] = {       4,       3.3,       2.4,       3.2,       4.3,         2};
   float cutpfhcalisooetom3[] = {    11.5,        22,       7.8,         8,       6.2,       7.4};
          float cutpfsieie3[] = {  0.0112,    0.0108,      0.01,     0.028,     0.029,     0.027};
         float cutpfhovere3[] = {   0.123,     0.142,      0.11,     0.146,     0.108,     0.059};
             float cutpfr93[] = {    0.94,       0.9,      0.24,      0.94,       0.9,      0.23};
  float cutpf_drtotk_25_993[] = {       1,         1,      0.06,    0.0163,     0.057,    0.0168};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.903277        fake=0.511431        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.000944444     end of iteration 8          HyperTight1  ---------------------------
    float cutpfisosumoet4[] = {     5.6,       4.8,       4.3,       5.6,       3.9,         3};
 float cutpfisosumoetbad4[] = {     9.8,       6.7,       5.6,       5.3,       4.6,         4};
       float cutpfchisooet4[] = {     3.4,       2.3,       2.4,       2.4,       1.9,      1.35};
   float cutpfhcalisooetom4[] = {    11.5,        22,       7.8,       4.3,       5.6,       3.3};
          float cutpfsieie4[] = {  0.0108,    0.0107,    0.0099,     0.028,     0.028,     0.027};
         float cutpfhovere4[] = {   0.123,     0.142,      0.11,     0.146,      0.09,     0.059};
             float cutpfr94[] = {    0.94,       0.9,      0.38,      0.94,       0.9,      0.23};
  float cutpf_drtotk_25_994[] = {       1,         1,      0.44,      0.36,      0.39,    0.0182};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.870999        fake=0.459021        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00122222     end of iteration 8          HyperTight2  ---------------------------
    float cutpfisosumoet5[] = {     5.2,       4.6,       3.6,         4,       3.2,       2.5};
 float cutpfisosumoetbad5[] = {       9,       6.4,       4.6,       4.9,       4.2,       3.8};
       float cutpfchisooet5[] = {     3.2,       2.1,       2.1,      1.97,      1.33,     0.028};
   float cutpfhcalisooetom5[] = {    11.5,        21,       7.8,       4.2,      1.93,       2.9};
          float cutpfsieie5[] = {  0.0107,    0.0107,    0.0098,     0.028,     0.027,     0.027};
         float cutpfhovere5[] = {   0.123,     0.142,      0.11,     0.108,     0.068,     0.059};
             float cutpfr95[] = {    0.94,       0.9,      0.38,      0.94,       0.9,      0.23};
  float cutpf_drtotk_25_995[] = {       1,         1,      0.99,      0.99,      0.96,     0.019};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.830285        fake=0.404951        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00166667     end of iteration 8          HyperTight3  ---------------------------
    float cutpfisosumoet6[] = {     4.8,       4.3,       2.7,       2.7,       2.5,       2.4};
 float cutpfisosumoetbad6[] = {     6.2,         6,       4.2,         4,       4.2,       3.4};
       float cutpfchisooet6[] = {     2.7,       2.1,      1.68,      1.25,       1.1,     0.028};
   float cutpfhcalisooetom6[] = {    11.5,        21,       7.6,       4.2,      1.93,      1.83};
          float cutpfsieie6[] = {  0.0107,    0.0105,    0.0098,     0.027,     0.027,     0.027};
         float cutpfhovere6[] = {   0.123,     0.122,      0.11,     0.108,     0.026,     0.059};
             float cutpfr96[] = {    0.94,       0.9,      0.38,      0.94,       0.9,      0.23};
  float cutpf_drtotk_25_996[] = {       1,         1,         1,         1,         1,      0.53};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.727129        fake=0.305297        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00222222     end of iteration 8          HyperTight4  ---------------------------
    float cutpfisosumoet7[] = {     2.8,       2.8,       2.4,       2.5,       2.4,       2.3};
 float cutpfisosumoetbad7[] = {     6.2,       5.4,       3.5,         4,       4.2,       2.6};
       float cutpfchisooet7[] = {    1.32,       1.3,      1.45,      0.76,      0.97,   0.00047};
   float cutpfhcalisooetom7[] = {     9.5,       4.3,       2.4,       4.2,      1.93,      1.34};
          float cutpfsieie7[] = {  0.0107,    0.0105,    0.0098,     0.027,     0.027,     0.027};
         float cutpfhovere7[] = {    0.09,     0.091,      0.11,     0.067,    0.0198,     0.053};
             float cutpfr97[] = {    0.94,       0.9,      0.38,      0.94,       0.9,      0.32};
  float cutpf_drtotk_25_997[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.534842        fake=0.174365        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00305556     end of iteration 8          HyperTight5  ---------------------------
    float cutpfisosumoet8[] = {     2.5,       2.4,       2.2,       2.3,       2.2,      1.98};
 float cutpfisosumoetbad8[] = {     5.7,       4.7,       2.7,         4,       4.2,       2.1};
       float cutpfchisooet8[] = {    1.07,         1,      1.34,      0.69,      0.69,     5e-06};
   float cutpfhcalisooetom8[] = {     9.5,       4.3,       2.4,       2.4,      1.93,      1.17};
          float cutpfsieie8[] = {  0.0107,    0.0104,    0.0097,     0.027,     0.024,     0.025};
         float cutpfhovere8[] = {   0.086,     0.091,     0.059,     0.057,    0.0198,     0.036};
             float cutpfr98[] = {    0.94,       0.9,      0.42,      0.95,       0.9,      0.53};
  float cutpf_drtotk_25_998[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.405524        fake=0.110157        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00416667     end of iteration 8          HyperTight6  ---------------------------
    float cutpfisosumoet9[] = {     2.4,       2.2,      1.88,       2.2,       2.2,       1.7};
 float cutpfisosumoetbad9[] = {     5.6,       3.6,       2.1,       2.4,       4.2,      1.55};
       float cutpfchisooet9[] = {    0.98,      0.68,      1.34,     0.007,      0.52,         0};
   float cutpfhcalisooetom9[] = {     5.6,       2.5,      1.27,         2,      1.02,      1.17};
          float cutpfsieie9[] = {  0.0107,    0.0103,    0.0097,     0.025,     0.023,     0.024};
         float cutpfhovere9[] = {   0.039,      0.04,     0.059,     0.032,    0.0198,    0.0178};
             float cutpfr99[] = {    0.94,       0.9,      0.42,      0.95,       0.9,      0.67};
  float cutpf_drtotk_25_999[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.288559        fake=0.0665108        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.00555556     end of iteration 8          HyperTight7  ---------------------------
   float cutpfisosumoet10[] = {     2.2,       2.1,      1.37,       2.1,       2.2,       1.7};
float cutpfisosumoetbad10[] = {     4.9,       3.6,      1.17,       2.1,       3.1,      1.55};
      float cutpfchisooet10[] = {    0.96,     0.149,    0.0135,     7e-05,     0.159,         0};
  float cutpfhcalisooetom10[] = {     2.1,      1.12,      0.69,      1.34,      0.78,       0.7};
         float cutpfsieie10[] = {  0.0107,    0.0101,     0.009,     0.024,     0.023,     0.023};
        float cutpfhovere10[] = {    0.02,    0.0164,     0.021,     0.029,     0.018,    0.0178};
            float cutpfr910[] = {    0.94,       0.9,      0.46,      0.95,       0.9,      0.76};
 float cutpf_drtotk_25_9910[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.193294        fake=0.0391879        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.0075     end of iteration 8          HyperTight8  ---------------------------
   float cutpfisosumoet11[] = {       2,      1.99,      1.35,         2,       2.2,       1.7};
float cutpfisosumoetbad11[] = {     3.8,       3.6,      1.03,       2.1,       3.1,      1.26};
      float cutpfchisooet11[] = {  0.0097,    0.0015,  0.000136,     1e-06,     0.159,         0};
  float cutpfhcalisooetom11[] = {    1.12,      0.23,      0.58,      1.14,    0.0078,       0.7};
         float cutpfsieie11[] = {  0.0099,      0.01,     0.009,     0.024,     0.023,     0.023};
        float cutpfhovere11[] = {   0.015,    0.0124,    0.0184,    0.0136,    0.0122,    0.0095};
            float cutpfr911[] = {    0.95,       0.9,      0.58,      0.95,       0.9,       0.8};
 float cutpf_drtotk_25_9911[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0573526        fake=0.0108523        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  ----------------------------sob value=0.01     end of iteration 8          HyperTight9  ---------------------------
   float cutpfisosumoet12[] = {       2,       1.8,      1.23,      1.89,       2.2,       1.7};
float cutpfisosumoetbad12[] = {     2.4,      1.65,      0.99,      1.58,       3.1,      1.26};
      float cutpfchisooet12[] = { 9.8e-05,   1.5e-05,     1e-06,         0,    0.0027,         0};
  float cutpfhcalisooetom12[] = {     0.6,      0.23,      0.58,      1.14,    0.0028,       0.7};
         float cutpfsieie12[] = {  0.0099,    0.0099,    0.0088,     0.024,     0.023,     0.023};
        float cutpfhovere12[] = {  0.0046,    0.0081,    0.0184,    0.0029,    0.0121,    0.0095};
            float cutpfr912[] = {    0.95,       0.9,      0.64,      0.95,      0.92,       0.8};
 float cutpf_drtotk_25_9912[] = {       1,         1,         1,         1,         1,         1};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0450285        fake=0.00877436        <<<<<<<<<<<<<<<<<<<<<<<<<<<<



  for(int ic=0; ic<6; ic++) {

    cutpfisosumoet6c[0][ic]=cutpfisosumoet0[ic];     cutpfisosumoet6c[1][ic]=cutpfisosumoet1[ic];     cutpfisosumoet6c[2][ic]=cutpfisosumoet2[ic];     cutpfisosumoet6c[3][ic]=cutpfisosumoet3[ic];    
    cutpfisosumoet6c[4][ic]=cutpfisosumoet4[ic];     cutpfisosumoet6c[5][ic]=cutpfisosumoet5[ic];     cutpfisosumoet6c[6][ic]=cutpfisosumoet6[ic];     cutpfisosumoet6c[7][ic]=cutpfisosumoet7[ic];
    cutpfisosumoet6c[8][ic]=cutpfisosumoet8[ic];     cutpfisosumoet6c[9][ic]=cutpfisosumoet9[ic];     cutpfisosumoet6c[10][ic]=cutpfisosumoet10[ic];   

    cutpfisosumoetbad6c[0][ic]=cutpfisosumoetbad0[ic];     cutpfisosumoetbad6c[1][ic]=cutpfisosumoetbad1[ic];     cutpfisosumoetbad6c[2][ic]=cutpfisosumoetbad2[ic];     cutpfisosumoetbad6c[3][ic]=cutpfisosumoetbad3[ic];    
    cutpfisosumoetbad6c[4][ic]=cutpfisosumoetbad4[ic];     cutpfisosumoetbad6c[5][ic]=cutpfisosumoetbad5[ic];     cutpfisosumoetbad6c[6][ic]=cutpfisosumoetbad6[ic];     cutpfisosumoetbad6c[7][ic]=cutpfisosumoetbad7[ic];
    cutpfisosumoetbad6c[8][ic]=cutpfisosumoetbad8[ic];     cutpfisosumoetbad6c[9][ic]=cutpfisosumoetbad9[ic];     cutpfisosumoetbad6c[10][ic]=cutpfisosumoetbad10[ic];   
    
    cutpfchisooet6c[0][ic]=cutpfchisooet0[ic];     cutpfchisooet6c[1][ic]=cutpfchisooet1[ic];     cutpfchisooet6c[2][ic]=cutpfchisooet2[ic];     cutpfchisooet6c[3][ic]=cutpfchisooet3[ic];    
    cutpfchisooet6c[4][ic]=cutpfchisooet4[ic];     cutpfchisooet6c[5][ic]=cutpfchisooet5[ic];     cutpfchisooet6c[6][ic]=cutpfchisooet6[ic];     cutpfchisooet6c[7][ic]=cutpfchisooet7[ic];
    cutpfchisooet6c[8][ic]=cutpfchisooet8[ic];     cutpfchisooet6c[9][ic]=cutpfchisooet9[ic];     cutpfchisooet6c[10][ic]=cutpfchisooet10[ic];   
    
    cutpfhcalisooetom6c[0][ic]=cutpfhcalisooetom0[ic];     cutpfhcalisooetom6c[1][ic]=cutpfhcalisooetom1[ic];     cutpfhcalisooetom6c[2][ic]=cutpfhcalisooetom2[ic];     cutpfhcalisooetom6c[3][ic]=cutpfhcalisooetom3[ic];    
    cutpfhcalisooetom6c[4][ic]=cutpfhcalisooetom4[ic];     cutpfhcalisooetom6c[5][ic]=cutpfhcalisooetom5[ic];     cutpfhcalisooetom6c[6][ic]=cutpfhcalisooetom6[ic];     cutpfhcalisooetom6c[7][ic]=cutpfhcalisooetom7[ic];
    cutpfhcalisooetom6c[8][ic]=cutpfhcalisooetom8[ic];     cutpfhcalisooetom6c[9][ic]=cutpfhcalisooetom9[ic];     cutpfhcalisooetom6c[10][ic]=cutpfhcalisooetom10[ic];   
    
    cutpfsieie6c[0][ic]=cutpfsieie0[ic];     cutpfsieie6c[1][ic]=cutpfsieie1[ic];     cutpfsieie6c[2][ic]=cutpfsieie2[ic];     cutpfsieie6c[3][ic]=cutpfsieie3[ic];    
    cutpfsieie6c[4][ic]=cutpfsieie4[ic];     cutpfsieie6c[5][ic]=cutpfsieie5[ic];     cutpfsieie6c[6][ic]=cutpfsieie6[ic];     cutpfsieie6c[7][ic]=cutpfsieie7[ic];
    cutpfsieie6c[8][ic]=cutpfsieie8[ic];     cutpfsieie6c[9][ic]=cutpfsieie9[ic];     cutpfsieie6c[10][ic]=cutpfsieie10[ic];   
    
    cutpfhovere6c[0][ic]=cutpfhovere0[ic];     cutpfhovere6c[1][ic]=cutpfhovere1[ic];     cutpfhovere6c[2][ic]=cutpfhovere2[ic];     cutpfhovere6c[3][ic]=cutpfhovere3[ic];    
    cutpfhovere6c[4][ic]=cutpfhovere4[ic];     cutpfhovere6c[5][ic]=cutpfhovere5[ic];     cutpfhovere6c[6][ic]=cutpfhovere6[ic];     cutpfhovere6c[7][ic]=cutpfhovere7[ic];
    cutpfhovere6c[8][ic]=cutpfhovere8[ic];     cutpfhovere6c[9][ic]=cutpfhovere9[ic];     cutpfhovere6c[10][ic]=cutpfhovere10[ic];   
    
    cutpfr96c[0][ic]=cutpfr90[ic];     cutpfr96c[1][ic]=cutpfr91[ic];     cutpfr96c[2][ic]=cutpfr92[ic];     cutpfr96c[3][ic]=cutpfr93[ic];    
    cutpfr96c[4][ic]=cutpfr94[ic];     cutpfr96c[5][ic]=cutpfr95[ic];     cutpfr96c[6][ic]=cutpfr96[ic];     cutpfr96c[7][ic]=cutpfr97[ic];
    cutpfr96c[8][ic]=cutpfr98[ic];     cutpfr96c[9][ic]=cutpfr99[ic];     cutpfr96c[10][ic]=cutpfr910[ic];   

    cutpf_drtotk_25_996c[0][ic]=cutpf_drtotk_25_990[ic];   cutpf_drtotk_25_996c[1][ic]=cutpf_drtotk_25_991[ic];   cutpf_drtotk_25_996c[2][ic]=cutpf_drtotk_25_992[ic];   cutpf_drtotk_25_996c[3][ic]=cutpf_drtotk_25_993[ic];    
    cutpf_drtotk_25_996c[4][ic]=cutpf_drtotk_25_994[ic];   cutpf_drtotk_25_996c[5][ic]=cutpf_drtotk_25_995[ic];   cutpf_drtotk_25_996c[6][ic]=cutpf_drtotk_25_996[ic];   cutpf_drtotk_25_996c[7][ic]=cutpf_drtotk_25_997[ic];
    cutpf_drtotk_25_996c[8][ic]=cutpf_drtotk_25_998[ic];   cutpf_drtotk_25_996c[9][ic]=cutpf_drtotk_25_999[ic];   cutpf_drtotk_25_996c[10][ic]=cutpf_drtotk_25_9910[ic];   

    cutpfisosumoet6c[11][ic]=cutpfisosumoet11[ic];     cutpfisosumoetbad6c[11][ic]=cutpfisosumoetbad11[ic];     cutpfchisooet6c[11][ic]=cutpfchisooet11[ic];     cutpfhcalisooetom6c[11][ic]=cutpfhcalisooetom11[ic];
    cutpfsieie6c[11][ic]=cutpfsieie11[ic];     cutpfhovere6c[11][ic]=cutpfhovere11[ic];     cutpfr96c[11][ic]=cutpfr911[ic];
    cutpf_drtotk_25_996c[11][ic]=cutpf_drtotk_25_9911[ic];   cutpf_drtotk_25_996c[11][ic]=cutpf_drtotk_25_9911[ic];

    cutpfisosumoet6c[12][ic]=cutpfisosumoet12[ic];     cutpfisosumoetbad6c[12][ic]=cutpfisosumoetbad12[ic];     cutpfchisooet6c[12][ic]=cutpfchisooet12[ic];     cutpfhcalisooetom6c[12][ic]=cutpfhcalisooetom12[ic];
    cutpfsieie6c[12][ic]=cutpfsieie12[ic];     cutpfhovere6c[12][ic]=cutpfhovere12[ic];     cutpfr96c[12][ic]=cutpfr912[ic];
    cutpf_drtotk_25_996c[12][ic]=cutpf_drtotk_25_9912[ic];   cutpf_drtotk_25_996c[12][ic]=cutpf_drtotk_25_9912[ic];

  }
}

void StatAnalysisExclusive::setPhotonCuts4() {
  //  Inits the cut levels in 4 categories for photon selection.  Just called once.
  //  Output of cut-setting program can just be copied here.

  //                           ----------------------------sob value=0.0002     end of iteration 8          Loose  ---------------------------
  float       cutsubleadisosumoet0[] = {      8.2,       4.1,       5.4,       2.6 };
  float    cutsubleadisosumoetbad0[] = {       67,        69,        85,       7.2 };
  float     cutsubleadtrkisooetom0[] = {      7.5,       4.5,       5.2,       2.5 };
  float           cutsubleadsieie0[] = {   0.0112,    0.0102,     0.029,     0.028 };
  float          cutsubleadhovere0[] = {     0.09,     0.089,     0.101,     0.073 };
  float              cutsubleadr90[] = {     0.94,      0.31,      0.92,      0.29 };
  float   cutsublead_drtotk_25_990[] = {     0.26,     0.029,    0.0062,    0.0055 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.940332        fake=0.0848119        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //                           ----------------------------sob value=0.0004     end of iteration 6          Medium  ---------------------------
  float       cutsubleadisosumoet1[] = {      6.4,       3.2,       3.4,       2.2 };
  float    cutsubleadisosumoetbad1[] = {       64,      10.8,        13,       3.5 };
  float     cutsubleadtrkisooetom1[] = {      6.4,       3.4,       3.8,       2.1 };
  float           cutsubleadsieie1[] = {   0.0109,      0.01,     0.029,     0.028 };
  float          cutsubleadhovere1[] = {    0.089,     0.079,      0.09,     0.061 };
  float              cutsubleadr91[] = {     0.94,      0.32,      0.94,      0.29 };
  float   cutsublead_drtotk_25_991[] = {     0.98,     0.029,    0.0109,    0.0111 };
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.918779        fake=0.0596949        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //			       ----------------------------sob value=0.0008     end of iteration 6          Tight  ---------------------------
  float       cutsubleadisosumoet2[] = {      4.7,       2.8,       2.5,      1.46 };
  float    cutsubleadisosumoetbad2[] = {       62,       5.2,       7.3,       2.5 };
  float     cutsubleadtrkisooetom2[] = {      4.7,       2.9,       3.8,      1.63 };
  float           cutsubleadsieie2[] = {   0.0107,    0.0099,     0.028,     0.027 };
  float          cutsubleadhovere2[] = {    0.087,     0.065,     0.087,      0.05 };
  float              cutsubleadr92[] = {     0.94,      0.34,      0.94,      0.29 };
  float   cutsublead_drtotk_25_992[] = {        1,     0.029,     0.021,     0.028 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.887763        fake=0.042122        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //			       ----------------------------sob value=0.0016     end of iteration 6          SuperTight  ---------------------------
  float       cutsubleadisosumoet3[] = {      3.8,       2.2,      1.77,      1.29 };
  float    cutsubleadisosumoetbad3[] = {     11.7,       3.4,       3.9,      1.84 };
  float     cutsubleadtrkisooetom3[] = {      3.5,       2.2,       2.3,      1.45 };
  float           cutsubleadsieie3[] = {   0.0106,    0.0097,     0.028,     0.027 };
  float          cutsubleadhovere3[] = {    0.082,     0.062,     0.065,     0.048 };
  float              cutsubleadr93[] = {     0.94,      0.36,      0.94,      0.32 };
  float   cutsublead_drtotk_25_993[] = {        1,     0.062,      0.97,      0.97 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.844291        fake=0.0290732        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //			       ----------------------------sob value=0.0032     end of iteration 6          HyperTight1  ---------------------------
  float       cutsubleadisosumoet4[] = {      3.2,      1.76,      1.39,      1.18 };
  float    cutsubleadisosumoetbad4[] = {      6.1,       2.7,       2.8,      0.66 };
  float     cutsubleadtrkisooetom4[] = {      3.4,      1.86,      1.67,      1.44 };
  float           cutsubleadsieie4[] = {   0.0104,    0.0094,     0.028,     0.025 };
  float          cutsubleadhovere4[] = {    0.076,      0.03,     0.047,     0.046 };
  float              cutsubleadr94[] = {     0.94,      0.41,      0.94,      0.34 };
  float   cutsublead_drtotk_25_994[] = {        1,      0.97,         1,         1 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.774013        fake=0.0190779        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //			       ----------------------------sob value=0.00625     end of iteration 6          HyperTight2  ---------------------------
  float       cutsubleadisosumoet5[] = {      2.6,      1.31,      1.33,      0.82 };
  float    cutsubleadisosumoetbad5[] = {      5.1,      1.62,      1.38,  -0.224864 };
  float     cutsubleadtrkisooetom5[] = {      2.9,       1.6,      1.55,      1.44 };
  float           cutsubleadsieie5[] = {   0.0101,    0.0093,     0.027,     0.023 };
  float          cutsubleadhovere5[] = {    0.048,    0.0189,     0.032,    0.0085 };
  float              cutsubleadr95[] = {     0.94,      0.47,      0.94,      0.52 };
  float   cutsublead_drtotk_25_995[] = {        1,         1,         1,         1 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.676391        fake=0.0121019        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //			       ----------------------------sob value=0.0125     end of iteration 6          HyperTight3  ---------------------------
  float       cutsubleadisosumoet6[] = {     1.85,      0.96,      1.21,  -0.028513 };
  float    cutsubleadisosumoetbad6[] = {      3.7,      0.97,      1.38,  -0.880416 };
  float     cutsubleadtrkisooetom6[] = {     1.93,       1.4,      1.48,     0.056 };
  float           cutsubleadsieie6[] = {   0.0099,    0.0092,     0.027,     0.023 };
  float          cutsubleadhovere6[] = {    0.042,    0.0173,     0.023,    0.0085 };
  float              cutsubleadr96[] = {     0.94,      0.69,      0.97,      0.52 };
  float   cutsublead_drtotk_25_996[] = {        1,         1,         1,         1 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.524271        fake=0.00631764        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //                            ----------------------------sob value=0.025     end of iteration 6          HyperTight4  ---------------------------
  float       cutsubleadisosumoet7[] = {     1.31,       0.3,      1.15,  -0.028513 };
  float    cutsubleadisosumoetbad7[] = {     1.72,      0.69,      1.14,  -0.880416 };
  float     cutsubleadtrkisooetom7[] = {     1.42,      0.76,      1.48,     0.056 };
  float           cutsubleadsieie7[] = {   0.0098,     0.009,     0.026,     0.023 };
  float          cutsubleadhovere7[] = {    0.037,   0.00049,    0.0198,   0.00024 };
  float              cutsubleadr97[] = {     0.94,      0.69,      0.97,      0.73 };
  float   cutsublead_drtotk_25_997[] = {        1,         1,         1,         1 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.383546        fake=0.00362626        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //                            ----------------------------sob value=0.045     end of iteration 6          HyperTight5  ---------------------------
  float       cutsubleadisosumoet8[] = {     0.94,     0.112,     0.039,  -0.028685 };
  float    cutsubleadisosumoetbad8[] = {     0.86,      0.45,  -0.371686,  -0.880416 };
  float     cutsubleadtrkisooetom8[] = {     1.21,      0.51,      0.27,     0.056 };
  float           cutsubleadsieie8[] = {   0.0097,     0.009,     0.026,     0.023 };
  float          cutsubleadhovere8[] = {    0.028,   1.4e-05,    0.0198,     7e-06 };
  float              cutsubleadr98[] = {     0.94,      0.69,      0.97,      0.73 };
  float   cutsublead_drtotk_25_998[] = {        1,         1,         1,         1 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.265878        fake=0.00218518        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //			         ----------------------------sob value=0.08     end of iteration 6          HyperTight6  ---------------------------
  float       cutsubleadisosumoet9[] = {      0.3,     0.072,  -0.039554,  -0.972925 };
  float    cutsubleadisosumoetbad9[] = {     0.59,      0.42,  -0.840197,  -1.76895 };
  float     cutsubleadtrkisooetom9[] = {     0.43,       0.4,    0.0075,     0.056 };
  float           cutsubleadsieie9[] = {   0.0094,     0.009,     0.024,     0.023 };
  float          cutsubleadhovere9[] = {   0.0071,         0,   0.00055,         0 };
  float              cutsubleadr99[] = {     0.95,      0.69,      0.97,      0.84 };
  float   cutsublead_drtotk_25_999[] = {        1,         1,         1,         1 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.109842        fake=0.000923024        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //			         ----------------------------sob value=0.15     end of iteration 6          HyperTight7  ---------------------------
  float       cutsubleadisosumoet10[] = {    0.146,     0.058,  -0.039554,  -1.16081 };
  float    cutsubleadisosumoetbad10[] = {     0.42,      0.41,  -0.840197,  -1.76895 };
  float     cutsubleadtrkisooetom10[] = {     0.43,    0.0111,    0.0075,    0.0073 };
  float           cutsubleadsieie10[] = {   0.0094,     0.009,     0.024,     0.023 };
  float          cutsubleadhovere10[] = {   0.0002,         0,   1.6e-05,         0 };
  float              cutsubleadr910[] = {     0.95,      0.69,      0.97,      0.85 };
  float   cutsublead_drtotk_25_9910[] = {        1,         1,         1,         1 };
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0898582        fake=0.00079285        <<<<<<<<<<<<<<<<<<<<<<<<<<<<

  for(int ic=0;ic!=4;++ic) {
    cutsubleadisosumoet[0][ic]=cutsubleadisosumoet0[ic];     cutsubleadisosumoet[1][ic]=cutsubleadisosumoet1[ic];     cutsubleadisosumoet[2][ic]=cutsubleadisosumoet2[ic];     cutsubleadisosumoet[3][ic]=cutsubleadisosumoet3[ic];    
    cutsubleadisosumoet[4][ic]=cutsubleadisosumoet4[ic];     cutsubleadisosumoet[5][ic]=cutsubleadisosumoet5[ic];     cutsubleadisosumoet[6][ic]=cutsubleadisosumoet6[ic];     cutsubleadisosumoet[7][ic]=cutsubleadisosumoet7[ic];
    cutsubleadisosumoet[8][ic]=cutsubleadisosumoet8[ic];     cutsubleadisosumoet[9][ic]=cutsubleadisosumoet9[ic];     cutsubleadisosumoet[10][ic]=cutsubleadisosumoet10[ic];   

    cutsubleadisosumoetbad[0][ic]=cutsubleadisosumoetbad0[ic];     cutsubleadisosumoetbad[1][ic]=cutsubleadisosumoetbad1[ic];     cutsubleadisosumoetbad[2][ic]=cutsubleadisosumoetbad2[ic];     cutsubleadisosumoetbad[3][ic]=cutsubleadisosumoetbad3[ic];    
    cutsubleadisosumoetbad[4][ic]=cutsubleadisosumoetbad4[ic];     cutsubleadisosumoetbad[5][ic]=cutsubleadisosumoetbad5[ic];     cutsubleadisosumoetbad[6][ic]=cutsubleadisosumoetbad6[ic];     cutsubleadisosumoetbad[7][ic]=cutsubleadisosumoetbad7[ic];
    cutsubleadisosumoetbad[8][ic]=cutsubleadisosumoetbad8[ic];     cutsubleadisosumoetbad[9][ic]=cutsubleadisosumoetbad9[ic];     cutsubleadisosumoetbad[10][ic]=cutsubleadisosumoetbad10[ic];   

    cutsubleadtrkisooetom[0][ic]=cutsubleadtrkisooetom0[ic];     cutsubleadtrkisooetom[1][ic]=cutsubleadtrkisooetom1[ic];     cutsubleadtrkisooetom[2][ic]=cutsubleadtrkisooetom2[ic];     cutsubleadtrkisooetom[3][ic]=cutsubleadtrkisooetom3[ic];    
    cutsubleadtrkisooetom[4][ic]=cutsubleadtrkisooetom4[ic];     cutsubleadtrkisooetom[5][ic]=cutsubleadtrkisooetom5[ic];     cutsubleadtrkisooetom[6][ic]=cutsubleadtrkisooetom6[ic];     cutsubleadtrkisooetom[7][ic]=cutsubleadtrkisooetom7[ic];
    cutsubleadtrkisooetom[8][ic]=cutsubleadtrkisooetom8[ic];     cutsubleadtrkisooetom[9][ic]=cutsubleadtrkisooetom9[ic];     cutsubleadtrkisooetom[10][ic]=cutsubleadtrkisooetom10[ic];   

    cutsubleadsieie[0][ic]=cutsubleadsieie0[ic];     cutsubleadsieie[1][ic]=cutsubleadsieie1[ic];     cutsubleadsieie[2][ic]=cutsubleadsieie2[ic];     cutsubleadsieie[3][ic]=cutsubleadsieie3[ic];    
    cutsubleadsieie[4][ic]=cutsubleadsieie4[ic];     cutsubleadsieie[5][ic]=cutsubleadsieie5[ic];     cutsubleadsieie[6][ic]=cutsubleadsieie6[ic];     cutsubleadsieie[7][ic]=cutsubleadsieie7[ic];
    cutsubleadsieie[8][ic]=cutsubleadsieie8[ic];     cutsubleadsieie[9][ic]=cutsubleadsieie9[ic];     cutsubleadsieie[10][ic]=cutsubleadsieie10[ic];   

    cutsubleadhovere[0][ic]=cutsubleadhovere0[ic];     cutsubleadhovere[1][ic]=cutsubleadhovere1[ic];     cutsubleadhovere[2][ic]=cutsubleadhovere2[ic];     cutsubleadhovere[3][ic]=cutsubleadhovere3[ic];    
    cutsubleadhovere[4][ic]=cutsubleadhovere4[ic];     cutsubleadhovere[5][ic]=cutsubleadhovere5[ic];     cutsubleadhovere[6][ic]=cutsubleadhovere6[ic];     cutsubleadhovere[7][ic]=cutsubleadhovere7[ic];
    cutsubleadhovere[8][ic]=cutsubleadhovere8[ic];     cutsubleadhovere[9][ic]=cutsubleadhovere9[ic];     cutsubleadhovere[10][ic]=cutsubleadhovere10[ic];   

    cutsubleadr9[0][ic]=cutsubleadr90[ic];     cutsubleadr9[1][ic]=cutsubleadr91[ic];     cutsubleadr9[2][ic]=cutsubleadr92[ic];     cutsubleadr9[3][ic]=cutsubleadr93[ic];    
    cutsubleadr9[4][ic]=cutsubleadr94[ic];     cutsubleadr9[5][ic]=cutsubleadr95[ic];     cutsubleadr9[6][ic]=cutsubleadr96[ic];     cutsubleadr9[7][ic]=cutsubleadr97[ic];
    cutsubleadr9[8][ic]=cutsubleadr98[ic];     cutsubleadr9[9][ic]=cutsubleadr99[ic];     cutsubleadr9[10][ic]=cutsubleadr910[ic];   

    cutsublead_drtotk_25_99[0][ic]=cutsublead_drtotk_25_990[ic];     cutsublead_drtotk_25_99[1][ic]=cutsublead_drtotk_25_991[ic];     cutsublead_drtotk_25_99[2][ic]=cutsublead_drtotk_25_992[ic];     cutsublead_drtotk_25_99[3][ic]=cutsublead_drtotk_25_993[ic];    
    cutsublead_drtotk_25_99[4][ic]=cutsublead_drtotk_25_994[ic];     cutsublead_drtotk_25_99[5][ic]=cutsublead_drtotk_25_995[ic];     cutsublead_drtotk_25_99[6][ic]=cutsublead_drtotk_25_996[ic];     cutsublead_drtotk_25_99[7][ic]=cutsublead_drtotk_25_997[ic];
    cutsublead_drtotk_25_99[8][ic]=cutsublead_drtotk_25_998[ic];     cutsublead_drtotk_25_99[9][ic]=cutsublead_drtotk_25_999[ic];     cutsublead_drtotk_25_99[10][ic]=cutsublead_drtotk_25_9910[ic];   

  }
  cout<<" setPhotonCuts4  "<<cutsubleadtrkisooetom[4][0]<<endl;				       
}

int StatAnalysisExclusive::photonCutLevel4(float r9, float eta, float et, float isosum, float isosumbad, float tkiso, float sieie, float hoe, float drtotk) {
  //  Returns the photon cut level for 4 categories from 0 (passing no cut level) to 11 (passing all cut levels).
  //  I think this should move to General Functions so that it can be used for photon ID in general, but also GA and CMSSW.

  if(MPDEBUGCOPY) cout << "photonCutLevel4 BEGIN"<<endl;
  if(MPDEBUGCOPY) cout<<" in photonCutLevel4  "<<r9<<"  "<<eta<<"  "<<isosum<<"  "<<isosumbad<<"  "<<tkiso<<"  "<<sieie<<"  "<<hoe<<"  "<<drtotk;

  if(!photonCutsSet4)   {setPhotonCuts4();     photonCutsSet4=true;}
  int ccat=0;
  if(fabs(eta)>1.479) ccat=2;                       //   4 cats
  if(r9<0.94) ccat++;
  int lev=0;
  while( lev < 11 && 
	 isosum*50/et    <= cutsubleadisosumoet[lev][ccat] &&
	 isosumbad*50/et <= cutsubleadisosumoetbad[lev][ccat] &&
	 tkiso*50/et     <= cutsubleadtrkisooetom[lev][ccat] &&
	 sieie           <= cutsubleadsieie[lev][ccat] &&
	 hoe             <= cutsubleadhovere[lev][ccat] &&
	 r9              >= cutsubleadr9[lev][ccat] &&
	 drtotk          >= cutsublead_drtotk_25_99[lev][ccat]                      )  lev++;

  if(MPDEBUGCOPY) cout<<"  "<<ccat<<"  "<<lev<<endl;
  if(MPDEBUGCOPY) cout << "PhotonCutLevel4 END"<<endl;
  return lev;

}

int StatAnalysisExclusive::photonCutLevel6(float r9, float eta, float et, float isosum, float isosumbad, float tkiso, float sieie, float hoe, float drtotk) {
  //  Returns the photon cut level for 6 categories from 0 (passing no cut level) to 11 (passing all cut levels).
  //  I think this should move to General Functions so that it can be used for photon ID in general, but also GA and CMSSW.

  if(MPDEBUGCOPY) cout << "photonCutLevel6 BEGIN"<<endl;

  if(!photonCutsSet6)   {setPhotonCuts6();     photonCutsSet6=true;}

  float etdepterm = 2.8;     // 5
  float etdeptermbad = 4.8;  // 7
  float etafac=0.5;
  float etafacbad=0.3;
  if(fabs(eta)>1.479) {etafac=0; etafacbad=0;}

  int ccat=0;
  if(fabs(eta)>1.479) ccat=3;                       //   4 cats
  if(r9<0.94) ccat++;
  if(r9<0.90) ccat++;
  int lev=0;

  while( lev < 11 && 
	 (isosum+etdepterm+etafac*fabs(eta))*50/et                  <= cutsubleadisosumoet6c[lev][ccat] &&
	 (isosumbad+etdeptermbad+etafacbad*fabs(eta))*50/et         <= cutsubleadisosumoetbad6c[lev][ccat] &&
	 tkiso*50/et                                                <= cutsubleadtrkisooetom6c[lev][ccat] &&
	 sieie                                                      <= cutsubleadsieie6c[lev][ccat] &&
	 hoe                                                        <= cutsubleadhovere6c[lev][ccat] &&
	 r9                                                         >= cutsubleadr96c[lev][ccat] &&
	 drtotk                                                     >= cutsublead_drtotk_25_996c[lev][ccat]                      )  lev++;

  if(MPDEBUGCOPY) cout << "PhotonCutLevel6 END"<<endl;
  return lev;

}

int StatAnalysisExclusive::photonCutLevel6pf(float r9, float eta, float et, float isosum, float isosumbad, float chargediso, float neutraliso, float sieie, float hoe, float drtotk) {
  //  Returns the photon cut level for 6 categories from 0 (passing no cut level) to 11 (passing all cut levels).
  //  I think this should move to General Functions so that it can be used for photon ID in general, but also GA and CMSSW.

  if(MPDEBUGCOPY) cout << "photonCutLevel6pf BEGIN"<<endl;

  if(!photonCutsSet6pf)   {setPhotonCuts6pf();     photonCutsSet6pf=true;}

  float etdepterm = 2.8;
  float etdeptermbad = 4.8;

  int ccat=0;
  if(fabs(eta)>1.479) ccat=3;                       //   4 cats
  if(r9<0.94) ccat++;
  if(r9<0.90) ccat++;
  int lev=0;

  while( lev < 13 && 
	 (isosum+etdepterm)*50/et               <= cutpfisosumoet6c[lev][ccat] &&
	 (isosumbad+etdeptermbad)*50/et         <= cutpfisosumoetbad6c[lev][ccat] &&
	 chargediso*50/et                       <= cutpfchisooet6c[lev][ccat] &&
	 neutraliso*50/et                       <= cutpfhcalisooetom6c[lev][ccat] &&
	 sieie                                  <= cutpfsieie6c[lev][ccat] &&
	 hoe                                    <= cutpfhovere6c[lev][ccat] &&
	 r9                                     >= cutpfr96c[lev][ccat] &&
	 drtotk                                 >= cutpf_drtotk_25_996c[lev][ccat]                      )  lev++;

  //cout<<"photonCutLevel6pf  lev="<<lev<<"  rho="<<rho<<"  r9="<<r9<<"  eta="<<eta<<"  et="<<et<<"  ccat="<<ccat<<"  isosum="<<isosum<<"  isosumbad="<<isosumbad<<"  chargediso="<<chargediso<<"  neutraliso="<<neutraliso
  //    <<"  sieie="<<sieie<<"  hoe="<<hoe<<"  drtotk="<<drtotk<<endl;
  //cout<<cutpfisosumoet6c[lev][ccat]<<"  "<<cutpfisosumoetbad6c[lev][ccat]<<"  "<<cutpfchisooet6c[lev][ccat]<<"  "<<cutpfhcalisooetom6c[lev][ccat]<<"  "<<cutpfsieie6c[lev][ccat]<<"  "<<cutpfhovere6c[lev][ccat]<<"  "
  //    <<cutpfr96c[lev][ccat]<<"  "<<cutpf_drtotk_25_996c[lev][ccat]<<endl;
  if(MPDEBUGCOPY) cout << "PhotonCutLevel6pf END"<<endl;
  return lev;
}


int StatAnalysisExclusive::diphoCutLevel(int leadind, int subleadind, int vtxind) {

  if(MPDEBUGCOPY) cout << "diphoCutLevel BEGIN"<<endl;
  //int leadind = dipho_leadind[diphoton_index];
  //int subleadind = dipho_subleadind[diphoton_index];

  int vtxi=vtxind, diphoind=-1;
  //  this loop doesn't work well because the two photons are not always in the list of diphotons
  if(vtxi<0) {
    //get right index to photons and variables....
    for(int idipho=0;idipho!=ll->dipho_n;++idipho) {
      if(MPDEBUG) cout<<"diphoCutLevel: leadind,subleadind"<<leadind<<" "<<subleadind<<" "<<ll->dipho_leadind[idipho]<<" "<<ll->dipho_subleadind[idipho]<<endl;
      if(leadind==ll->dipho_leadind[idipho]&&subleadind==ll->dipho_subleadind[idipho]) {
	vtxi=ll->dipho_vtxind[idipho];
	diphoind=idipho;
      }
      //     jgb-kludge for photon pairs that don't appear in the diphoton list
      //else if(diphoind==-1 && dipho_vtxind[idipho]==vtxind) 	diphoind=idipho;                         //  use some reasonable vertex if this diphoton does not appear in the list
    }
  }




//MARCO FIX

  //TLorentzVector newleadp4; //MARCO FIX=  TLorentzVector(VertexCorrectedP4Hgg(leadind,vtxind));
  //TLorentzVector newsubleadp4; //MARCO FIX=  TLorentzVector(VertexCorrectedP4Hgg(subleadind,vtxind));
//MARCO FIXED
  TLorentzVector newleadp4 = ll->get_pho_p4(leadind,vtxind);
  TLorentzVector newsubleadp4 = ll->get_pho_p4(subleadind,vtxind);

  TLorentzVector * leadp4 = &newleadp4;
  TLorentzVector * subleadp4 = &newsubleadp4;



  TVector3 leadcalopos = *((TVector3*)ll->pho_calopos->At(leadind));
  Float_t leadeta = fabs(((TLorentzVector*)ll->sc_p4->At(ll->pho_scind[leadind]))->Eta());

  TVector3 subleadcalopos = *((TVector3*)ll->pho_calopos->At(subleadind));
  Float_t subleadeta = fabs(((TLorentzVector*)ll->sc_p4->At(ll->pho_scind[subleadind]))->Eta());

  TLorentzVector diphotonp4 = (*leadp4) + (*subleadp4);

  //cout<<"cutdiphoptom="<<endl;
  //for(int j=0; j<4; j++) {for(int i=0; i<13; i++) {cout<<cutdiphoptom[i][j]<<"  ";} cout<<endl;}

  // list of variables needed below
  // diphopt
  // mass
  // (sub)leadcutindex
  // etamax
  // t_etamin
  // (sub)leadr9
  // (sub)t_leadpt
  float diphopt = diphotonp4.Pt();
  float mass = diphotonp4.M();


  //sean: this has to be changed too
  //MARCO FIXED OK???
  Float_t leadcutindex= PhotonCiCSelectionLevelJim(leadind,6,vtxi,0,diphoind);
  Float_t subleadcutindex= PhotonCiCSelectionLevelJim(subleadind,6,vtxi,1,diphoind);

  //cout<<"mmmmm "<<leadcutindex<<" "<<subleadcutindex<<endl;

  //std::pair<int,int> diphoton_indices(DiphotonCiCSelectionIndices( LEADCUTLEVEL, SUBLEADCUTLEVEL, leadPtMin, subleadPtMin, CICNCAT, -1, false));
  //leadind = diphoton_indices.first;
  //subleadind = diphoton_indices.second;


  float etamax = TMath::Max(leadeta,subleadeta);
  float etamin = TMath::Min(leadeta,subleadeta);
  float leadr9 = ll->pho_r9[leadind];
  float subleadr9 = ll->pho_r9[subleadind];
  float leadpt = leadp4->Pt();
  float subleadpt = subleadp4->Pt();

  float diphoptom=diphopt/mass;
  float fsubleadcutindex=float(11-subleadcutindex);
  float fleadcutindex=float(11-leadcutindex);
  //float subleadabseta=fabs(subleadeta);
  //float leadabseta=fabs(leadeta);
  //  float etamin=fabs(subleadeta); 
  //  float etamax-fabs(leadeta);
  //  if(etamin>etamax) {etamin=etamax;  etamax=fabs(subleadeta);}
  int dcat=0;
  if(etamax>1.479) dcat=2;
  if(subleadr9<0.94||leadr9<0.94) dcat++;
  float leadptomass=leadpt/mass;  
  float subleadptomass=subleadpt/mass; 
  float sumptom=subleadptomass+leadptomass; 

  float dmom=0.5*pow(pow(ll->pho_regr_energyerr[leadind]/ll->pho_regr_energy[leadind],2)+pow(ll->pho_regr_energyerr[subleadind]/ll->pho_regr_energy[subleadind],2),0.5);

  //cout<<"diphoCutLevel  "<<subleadind<<leadind<<vtxind<<vtxi<<diphoind<<"  "<<dcat<<"  "<<diphoptom<<"  "<<fsubleadcutindex<<"  "<<fleadcutindex<<"  "<<etamax<<"  "<<etamin<<"  "<<subleadptomass<<"  "<<leadptomass<<endl;

  int lev=0;
  //cout<<(bool)(diphoptom>cutdiphoptom[lev][dcat])<<(bool)(fsubleadcutindex<cutfsubleadcutindex[lev][dcat])<<(bool)(fleadcutindex<cutfleadcutindex[lev][dcat])<<(bool)(etamax<cutetamax[lev][dcat])
  //    <<(bool)(etamin<cutetamin[lev][dcat])<<(bool)(subleadptomass>cutsubleadptomass[lev][dcat])<<(bool)(leadptomass>cutleadptomass[lev][dcat])<<(bool)(mtom>cutmtom[lev][dcat])<<endl;

  while( 
	lev              < 14 && 
	diphoptom        > cutdiphoptom[lev][dcat] &&
	dmom             < cutdmom[lev][dcat] &&
	fsubleadcutindex < cutfsubleadcutindex[lev][dcat] &&
	fleadcutindex    < cutfleadcutindex[lev][dcat] &&
	etamax           < cutetamax[lev][dcat] &&
	subleadptomass   > cutsubleadptomass[lev][dcat] &&
	sumptom          > cutsumptom[lev][dcat])                    lev++;

  if(MPDEBUGCOPY) cout << "diphoCutLevel END"<<endl;
  return lev;

}

int StatAnalysisExclusive::diphoSubCategory(int nSubCat, int c4, int diphoCutLev) {
  if(nSubCat<=0) return 1;                                          //  no diphoton cuts if nSubCat<=0
  if(nSubCat>4) { cout<<"  nSubCat>4"<<endl;  return 1;}
  
  
  //     for 4 subcategories                                         disc      excl        cat0    cat1    cat2   cat3            diphoptom                           subleadcut
  //                                                        llrSB=-35.94     llrB=29.27   (17.48, -13.78, -2.66, -2.23)     0      1      2      3
  //                                                 data   llrSB=-32.46     llrB=27.42   (14.59   13.50   2.47   2.01)
  //int ic1a[]={ 3,  3,  1, 0};         //                                                                                0.0027  0.0066   0.038   0.026          T     ST      T      T
  //int ic2a[]={ 9,  5,  4, 2};         //                                                                                0.076   0.123    0.074   0.104        HT2    HT3    HT2    HT1
  //int ic3a[]={11,  8,  6, 5};         //                                                                                0.076   0.123    0.074   0.104        HT2    HT3    HT2    HT1
  //int ic4a[]={13, 13, 10, 6};         //                                                                                0.61    0.41     0.109   0.167        HT2    HT3    HT2    HT1

  //int ic1a[]={ 3,  1,  1, 0};         //                                                                                0.0027  0.0066   0.038   0.026          T     ST      T      T
  //int ic2a[]={ 5,  7,  4, 3};         //                                                                                0.076   0.123    0.074   0.104        HT2    HT3    HT2    HT1
  //int ic3a[]={ 8,  9,  6, 5};         //                                                                                0.076   0.123    0.074   0.104        HT2    HT3    HT2    HT1
  //int ic4a[]={13, 13,  9, 6};         //                                                                                0.61    0.41     0.109   0.167        HT2    HT3    HT2    HT1
  int ic0a[]={ 0,  0,  0, 0};     
  int ic1a[]={ 3,  1,  1, 0};         //                                                                                0.0027  0.0066   0.038   0.026          T     ST      T      T
  int ic2a[]={ 5,  7,  4, 3};         //                                                                                0.076   0.123    0.074   0.104        HT2    HT3    HT2    HT1
  int ic3a[]={ 8,  9,  6, 5};         //                                                                                0.076   0.123    0.074   0.104        HT2    HT3    HT2    HT1
  int ic4a[]={13, 13,  9, 6};         //                                                                                0.61    0.41     0.109   0.167        HT2    HT3    HT2    HT1
  int ic5a[]={15, 15, 15, 15};  
  //    for 3 subcategories                                         disc      excl        cat0    cat1    cat2   cat3            diphoptom                           subleadcut
  //                                                      mc       -38.17    30.21      (-18.73, -14.75, -2.61, -2.19)     0      1       2       3             0      1      2      3
  //                                                      data     -32.19    27.20      (-14.53, -13.41, -2.39, -2.00)
  //                                                               -23.2     19.41       ------  ------
  if(nSubCat==3) {
    //ic1a[0]=3;   ic1a[1]=1;   ic1a[2]=1;   ic1a[3]=0;         
    //ic2a[0]=8;   ic2a[1]=8;   ic2a[2]=5;   ic2a[3]=3;         
    //ic3a[0]=13;  ic3a[1]=13;  ic3a[2]=10;  ic3a[3]=6;         
    //ic4a[0]=15;  ic4a[1]=15;  ic4a[2]=15;  ic4a[3]=15;
    ic1a[0]=3;   ic1a[1]=1;   ic1a[2]=1;   ic1a[3]=0;         
    ic2a[0]=8;   ic2a[1]=8;   ic2a[2]=5;   ic2a[3]=3;         
    ic3a[0]=13;  ic3a[1]=13;  ic3a[2]=10;  ic3a[3]=6;         
    ic4a[0]=15;  ic4a[1]=15;  ic4a[2]=15;  ic4a[3]=15;
  }
  if(nSubCat==2) {
    //ic1a[0]=3;   ic1a[1]=3;   ic1a[2]=1;   ic1a[3]=1;           
    //ic2a[0]=10;   ic2a[1]=11;   ic2a[2]=5;   ic2a[3]=6;          -25.12    19.83     (-12.05, -10.21, -1.45,  -1.48)
    //                                                       data  -30.66    26.04     (-14.30, -12.33, -2.26,  -1.88)
    ic1a[0]=3;   ic1a[1]=3;   ic1a[2]=1;   ic1a[3]=1;         //                                        ------  -----   0.0027  0.037    0.038   0.026          T     ST      T      T
    ic2a[0]=10;  ic2a[1]=11;  ic2a[2]=5;   ic2a[3]=6;         //                                                        0.61    0.41     0.109   0.167        HT2    HT3    HT2    HT1
    ic3a[0]=15;  ic3a[1]=15;  ic3a[2]=15;  ic3a[3]=15;
  }
  if(nSubCat==1) {
    //ic1a[0]=3;   ic1a[1]=4;   ic1a[2]=2;   ic1a[3]=2;            -19.31    16.73    (-9.573,  -7.38, -1.209, -1.219)
    //                                         data             -27.92     llrB=24.2  (-13.21, -11.23, -1.92,  -1.66)
    ic1a[0]=3;   ic1a[1]=4;   ic1a[2]=2;   ic1a[3]=2;         //                                                        0.05    0.037    0.038   0.026        HT1     ST      T      T
    ic2a[0]=15;  ic2a[1]=15;  ic2a[2]=15;  ic2a[3]=15;
  }
  //   for no dipho selection                      mc      llrSB=-28.93     llrB=24.6
  //                                               data    llrSB=-27.28     llrB=23.5 (-12.57, -11.23, -1.92, -1.66)
  int *ic1[5];
  ic1[0]=ic0a;   ic1[1]=ic1a;   ic1[2]=ic2a;   ic1[3]=ic3a;   ic1[4]=ic4a;     ic1[5]=ic5a;  
  
  int diphoSubcat=0;
  for(int ic=0; ic<nSubCat+1; ic++) {if(diphoCutLev>=ic1[ic][c4]&&diphoCutLev<ic1[ic+1][c4]){diphoSubcat=ic; break;}}
  return diphoSubcat;
}

float StatAnalysisExclusive::diphoSubCategory(int c4, int diphoCutLev) {
  //  Returns the weight (s/b) for the diphoton.  Used to make weighted (mass) histograms in which the categories can be combined
  //  This weight includes the expected resolution in the category c4.
  //  Jim should update this as we get more data since s/b is measured from the data.
  float owsob[4][15]={{0,          0.393897,  0.433008,  0.433008,  0.836777,  1.26343,   1.24698,   1.21773,   1.66486,  1.58613,  1.67371,  2.05731,  1.92396,  2.70931,  3.19341},
		      {0,          0.245151,  0.301681,  0.468778,  0.608091,  0.654541,  0.824609,  0.955908,  1.14663,  1.14663,  1.14663,  1.49005,  1.08435,  2.15387,  2.06251},
		      {0,          0.103445,  0.180206,  0.214078,  0.308207,  0.396138,  0.515793,  0.509308,  0.403491, 1.11622,  1.11622,  1.11622,  1.11622,  1.11622,  1.11622},
		      {0.0853153,  0.111469,  0.18104,   0.239611,  0.297499,  0.33897,   0.518697,  0.518697,  0.518697, 0.518697, 0.518697, 0.518697, 0.518697, 0.518697, 0.518697}};
  return owsob[c4][diphoCutLev];
}

// Jim's Functions END




//BDT WAS IN GEN FUNCTIONS

void StatAnalysisExclusive::SetBDT() {

  //MPDEBUG=1;

  /*
  std::cout<<"SetBDT"<<std::endl;
  if(MPDEBUG)  std::cout<<"SetBDT"<<std::endl;

  tmvaReader = new TMVA::Reader("!Color:Silent");
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("sieie", &tmva_sieie);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("badpf_iso", &tmva_badpf_iso);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("drtotk", &tmva_drtotk);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("hoe", &tmva_hoe);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("tkisopf", &tmva_tkisopf);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("r9", &tmva_r9);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("pt", &tmva_pt);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddVariable("eta", &tmva_eta);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->AddSpectator("isLeading", &tmva_isLeading);
  std::cout<<"SetBDT"<<std::endl;
  tmvaReader->BookMVA("Gradient", "/home/users/branson/globeArea11/globe/mh_110_135_Gradient.weights.xml");
  std::cout<<"SetBDT"<<std::endl;


  tmvaReader1 = new TMVA::Reader("!Color:Silent");
  tmvaReader1->AddVariable("sieie", &tmva_sieie);
  tmvaReader1->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  tmvaReader1->AddVariable("badpf_iso", &tmva_badpf_iso);
  tmvaReader1->AddVariable("drtotk", &tmva_drtotk);
  tmvaReader1->AddVariable("hoe", &tmva_hoe);
  tmvaReader1->AddVariable("tkisopf", &tmva_tkisopf);
  tmvaReader1->AddVariable("r9", &tmva_r9);
  tmvaReader1->AddVariable("pt", &tmva_pt);
  tmvaReader1->AddVariable("eta", &tmva_eta);
  tmvaReader1->AddSpectator("isLeading", &tmva_isLeading);
  tmvaReader1->BookMVA("Category_Gradient", "/home/users/branson/globeArea11/globe/MVAweights/mh_110_135_Category_Gradient.weights.xml");

  std::cout<<"SetBDT1"<<std::endl;
  */
  tmvaReader2 = new TMVA::Reader("!Color:Silent");
  tmvaReader2->AddVariable("sieie", &tmva_sieie);
  tmvaReader2->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  tmvaReader2->AddVariable("badpf_iso", &tmva_badpf_iso);
  tmvaReader2->AddVariable("drtotk", &tmva_drtotk);
  tmvaReader2->AddVariable("hoe", &tmva_hoe);
  tmvaReader2->AddVariable("tkisopf", &tmva_tkisopf);
  tmvaReader2->AddVariable("r9", &tmva_r9);
  tmvaReader2->AddVariable("ptom", &tmva_ptom);
  tmvaReader2->AddVariable("eta", &tmva_eta);
  tmvaReader2->AddSpectator("isLeading", &tmva_isLeading);
  tmvaReader2->BookMVA("Gradient", "../ID_UCSD.weights.xml");
  if(MPDEBUG)  std::cout<<"SetBDT End"<<std::endl;

  /*
  tmvaReader3 = new TMVA::Reader("!Color:Silent");
  tmvaReader3->AddVariable("sieie", &tmva_sieie);
  tmvaReader3->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  tmvaReader3->AddVariable("badpf_iso", &tmva_badpf_iso);
  tmvaReader3->AddVariable("drtotk", &tmva_drtotk);
  tmvaReader3->AddVariable("hoe", &tmva_hoe);
  tmvaReader3->AddVariable("tkisopf", &tmva_tkisopf);
  tmvaReader3->AddVariable("r9", &tmva_r9);
  tmvaReader3->AddVariable("ptom", &tmva_ptom);
  tmvaReader3->AddVariable("eta", &tmva_eta);
  tmvaReader3->AddSpectator("isLeading", &tmva_isLeading);
  tmvaReader3->BookMVA("Gradient", "/home/users/matteo/test_Gradient.weights.xml");
  if(MPDEBUG)  std::cout<<"SetBDT End"<<std::endl;
  */

  tmvaReader_dipho = new TMVA::Reader("!Color:Silent"); 
  tmvaReader_dipho->AddVariable("subleadptomass", &tmva_subleadptomass);
  tmvaReader_dipho->AddVariable("diphoptom", &tmva_diphoptom);
  tmvaReader_dipho->AddVariable("sumptom", &tmva_sumptom);
  tmvaReader_dipho->AddVariable("subleadmva", &tmva_subleadmva);
  tmvaReader_dipho->AddVariable("leadmva", &tmva_leadmva);
  tmvaReader_dipho->AddVariable("leadeta", &tmva_leadeta);
  tmvaReader_dipho->AddVariable("subleadeta", &tmva_subleadeta);
  tmvaReader_dipho->AddVariable("leadr9", &tmva_leadr9);
  tmvaReader_dipho->AddVariable("subleadr9", &tmva_subleadr9);
  tmvaReader_dipho->AddVariable("dmom", &tmva_dmom);
  tmvaReader_dipho->AddSpectator("diphocat2r92eta", &tmva_isLeading);
  tmvaReader_dipho->BookMVA("Gradient", "../diphoton_UCSD.weights.xml");

  /*
  tmvaReader_dipho2 = new TMVA::Reader("!Color:Silent"); 
  tmvaReader_dipho2->AddVariable("subleadptomass", &tmva_subleadptomass);
  tmvaReader_dipho2->AddVariable("diphoptom", &tmva_diphoptom);
  tmvaReader_dipho2->AddVariable("sumptom", &tmva_sumptom);
  tmvaReader_dipho2->AddVariable("subleadmva", &tmva_subleadmva);
  tmvaReader_dipho2->AddVariable("leadmva", &tmva_leadmva);
  tmvaReader_dipho2->AddVariable("leadeta", &tmva_leadeta);
  tmvaReader_dipho2->AddVariable("subleadeta", &tmva_subleadeta);
  tmvaReader_dipho2->AddVariable("leadr9", &tmva_leadr9);
  tmvaReader_dipho2->AddVariable("subleadr9", &tmva_subleadr9);
  tmvaReader_dipho2->AddSpectator("diphocat2r92eta", &tmva_isLeading);
  tmvaReader_dipho2->BookMVA("Gradient", "/home/users/matteo/diphoton_allstat_Gradient.weights.xml");
  */
}


Float_t StatAnalysisExclusive::BDT(Int_t jentry, Int_t iPhoton, Int_t vtx) {
 
  if(MPDEBUG)  std::cout<<"BDT"<<std::endl;

  int n_r9_categories = 3;
  int n_eta_categories = 2;
  tmva_cat = ll->PhotonCategory(iPhoton,n_r9_categories,n_eta_categories);

  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;
  ll->BdtGetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 after"<<std::endl;

  /*
  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;

  if(jentry>0) {
  if (b_pho_sieie->GetReadEntry() != jentry)
  b_pho_sieie->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 1"<<std::endl;
  if (b_pho_pfiso_mycharged04->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 2"<<std::endl;
  if (b_pho_pfiso_myphoton04->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 3"<<std::endl;
  if (b_pho_pfiso_mycharged03->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged03->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 4"<<std::endl;
  if (b_pho_pfiso_myphoton03->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton03->GetEntry(jentry);
  //std::cout<<"BDT 1 - 5"<<std::endl;
  //if (b_pho_drtotk_25_99->GetReadEntry() != jentry)
  //b_pho_drtotk_25_99->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 6"<<std::endl;
  if (b_pho_hoe->GetReadEntry() != jentry)
  b_pho_hoe->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 7"<<std::endl;
  if (b_pho_r9->GetReadEntry() != jentry)
  b_pho_r9->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 8"<<std::endl;
  if (b_pho_p4->GetReadEntry() != jentry)
  b_pho_p4->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 9"<<std::endl;
  if (b_rho->GetReadEntry() != jentry)
  b_rho->GetEntry(jentry);
  }

  if(MPDEBUG)  std::cout<<"BDT 2"<<std::endl;
  */



  float isomax=-99;   int badind=0;
  for(int iv=0; iv<ll->vtx_std_n; iv++) if((*ll->pho_pfiso_mycharged04)[iPhoton][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged04)[iPhoton][iv]; }
  if(MPDEBUG)  std::cout<<"BDT 3"<<std::endl;
  
  float rhofacpf[6]={0.075, 0.082, 0.143, 0.050, 0.091, 0.106};          //move
  float rhofacbadpf[6]={0.141, 0.149, 0.208, 0.135, 0.162, 0.165};
  float rhofac=rhofacpf[tmva_cat];
  float rhofacbad=rhofacbadpf[tmva_cat];

  tmva_pt = ((TLorentzVector*)ll->pho_p4->At(iPhoton))->Et();
  //tmva_badpf_iso = ((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad+4.8)*50/tmva_pt;
  //tmva_goodpf_iso = ((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac+2.8)*50/tmva_pt;
  tmva_badpf_iso = ((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad)*50/tmva_pt;
  tmva_goodpf_iso = ((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac)*50/tmva_pt;
  tmva_tkisopf = (*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt;
  if(MPDEBUG)  std::cout<<"BDT 4"<<std::endl;
  
  tmva_sieie = ll->pho_sieie[iPhoton];
  tmva_drtotk = ll->pho_drtotk_25_99[iPhoton];
  tmva_hoe = ll->pho_hoe[iPhoton];
  tmva_r9 = ll->pho_r9[iPhoton];
  tmva_eta = fabs(((TLorentzVector*)ll->pho_p4->At(iPhoton))->Eta());
  
  if(MPDEBUG)  std::cout<<"BDT end"<<std::endl;
  
  Float_t mva = tmvaReader->EvaluateMVA("Gradient");
  //cout<<"BDT  "<<jentry<<"  "<<iPhoton<<"  "<<vtx<<"  "<<isomax<<"  "<<badind<<"  "<<rhofac<<"  "<<rhofacbad<<"  "<<tmva_pt<<"  "<<tmva_badpf_iso<<"  "<<tmva_goodpf_iso<<"  "<<tmva_sieie<<"  "<<tmva_drtotk
  //     <<"  "<<tmva_hoe<<"  "<<tmva_tkisopf<<"  "<<tmva_r9<<"  "<<tmva_eta<<"  "<<mva<<endl;
  return mva;
}

Float_t StatAnalysisExclusive::BDT_categorized(Int_t jentry, Int_t iPhoton, Int_t vtx, float isLead) {
 
  if(MPDEBUG)  std::cout<<"BDT"<<std::endl;

  int n_r9_categories = 3;
  int n_eta_categories = 2;
  tmva_cat = ll->PhotonCategory(iPhoton,n_r9_categories,n_eta_categories);

  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;
  ll->BdtGetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 after"<<std::endl;
  /*
  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;

  if(jentry>0) {
  if (b_pho_sieie->GetReadEntry() != jentry)
  b_pho_sieie->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 1"<<std::endl;
  if (b_pho_pfiso_mycharged04->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 2"<<std::endl;
  if (b_pho_pfiso_myphoton04->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 3"<<std::endl;
  if (b_pho_pfiso_mycharged03->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged03->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 4"<<std::endl;
  if (b_pho_pfiso_myphoton03->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton03->GetEntry(jentry);
  //std::cout<<"BDT 1 - 5"<<std::endl;
  //if (b_pho_drtotk_25_99->GetReadEntry() != jentry)
  //b_pho_drtotk_25_99->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 6"<<std::endl;
  if (b_pho_hoe->GetReadEntry() != jentry)
  b_pho_hoe->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 7"<<std::endl;
  if (b_pho_r9->GetReadEntry() != jentry)
  b_pho_r9->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 8"<<std::endl;
  if (b_pho_p4->GetReadEntry() != jentry)
  b_pho_p4->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 9"<<std::endl;
  if (b_rho->GetReadEntry() != jentry)
  b_rho->GetEntry(jentry);
  }

  if(MPDEBUG)  std::cout<<"BDT 2"<<std::endl;
  */


  float isomax=-99;   int badind=0;
  for(int iv=0; iv<ll->vtx_std_n; iv++) if((*ll->pho_pfiso_mycharged04)[iPhoton][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged04)[iPhoton][iv]; }
  if(MPDEBUG)  std::cout<<"BDT 3"<<std::endl;
  
  float rhofacpf[6]={0.075, 0.082, 0.143, 0.050, 0.091, 0.106};          //move
  float rhofacbadpf[6]={0.141, 0.149, 0.208, 0.135, 0.162, 0.165};
  float rhofac=rhofacpf[tmva_cat];
  float rhofacbad=rhofacbadpf[tmva_cat];

  tmva_pt = ((TLorentzVector*)ll->pho_p4->At(iPhoton))->Et();
  tmva_isLeading = isLead;
  //tmva_badpf_iso = ((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad+4.8)*50/tmva_pt;
  //tmva_goodpf_iso = ((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac+2.8)*50/tmva_pt;
  tmva_badpf_iso = ((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad)*50/tmva_pt;
  tmva_goodpf_iso = ((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac)*50/tmva_pt;
  tmva_tkisopf = (*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt;
  if(MPDEBUG)  std::cout<<"BDT 4"<<std::endl;
  
  tmva_sieie = ll->pho_sieie[iPhoton];
  tmva_drtotk = ll->pho_drtotk_25_99[iPhoton];
  tmva_hoe = ll->pho_hoe[iPhoton];
  tmva_r9 = ll->pho_r9[iPhoton];
  tmva_eta = fabs(((TLorentzVector*)ll->pho_p4->At(iPhoton))->Eta());
  
  if(MPDEBUG)  std::cout<<"BDT end"<<std::endl;
  
  Float_t mva = tmvaReader1->EvaluateMVA("Category_Gradient");
  //cout<<"BDT  "<<jentry<<"  "<<iPhoton<<"  "<<vtx<<"  "<<isomax<<"  "<<badind<<"  "<<rhofac<<"  "<<rhofacbad<<"  "<<tmva_pt<<"  "<<tmva_badpf_iso<<"  "<<tmva_goodpf_iso<<"  "<<tmva_sieie<<"  "<<tmva_drtotk
  //    <<"  "<<tmva_hoe<<"  "<<tmva_tkisopf<<"  "<<tmva_r9<<"  "<<tmva_eta<<"  "<<mva<<endl;
  return mva;
}

Float_t StatAnalysisExclusive::BDT_ptom(Int_t jentry, Int_t iPhoton, Int_t vtx, float mass) {

  // we put in the diphoton mass just to scale the single photon pt (pt/m) 
  if(MPDEBUG)  std::cout<<"BDT"<<std::endl;

  int n_r9_categories = 3;
  int n_eta_categories = 2;
  tmva_cat = ll->PhotonCategory(iPhoton,n_r9_categories,n_eta_categories);

  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;
  ll->BdtGetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 after"<<std::endl;

  /*
  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;

  if (b_pho_sieie->GetReadEntry() != jentry)
  b_pho_sieie->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 1"<<std::endl;
  if (b_pho_pfiso_mycharged04->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 2"<<std::endl;
  if (b_pho_pfiso_myphoton04->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 3"<<std::endl;
  if (b_pho_pfiso_mycharged03->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged03->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 4"<<std::endl;
  if (b_pho_pfiso_myphoton03->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton03->GetEntry(jentry);
  //std::cout<<"BDT 1 - 5"<<std::endl;
  //if (b_pho_drtotk_25_99->GetReadEntry() != jentry)
  //b_pho_drtotk_25_99->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 6"<<std::endl;
  if (b_pho_hoe->GetReadEntry() != jentry)
  b_pho_hoe->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 7"<<std::endl;
  if (b_pho_r9->GetReadEntry() != jentry)
  b_pho_r9->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 8"<<std::endl;
  if (b_pho_p4->GetReadEntry() != jentry)
  b_pho_p4->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 9"<<std::endl;
  if (b_rho->GetReadEntry() != jentry)
  b_rho->GetEntry(jentry);
  */

  if(MPDEBUG)  std::cout<<"BDT 2"<<std::endl;
  
  float isomax=-99;   int badind=0;
  for(int iv=0; iv<ll->vtx_std_n; iv++) if((*ll->pho_pfiso_mycharged04)[iPhoton][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged04)[iPhoton][iv]; }
  if(MPDEBUG)  std::cout<<"BDT 3"<<std::endl;
  
  float rhofacpf[6]={0.075, 0.082, 0.143, 0.050, 0.091, 0.106};          //move
  float rhofacbadpf[6]={0.141, 0.149, 0.208, 0.135, 0.162, 0.165};
  float rhofac=rhofacpf[tmva_cat];
  float rhofacbad=rhofacbadpf[tmva_cat];

  tmva_pt = ((TLorentzVector*)ll->pho_p4->At(iPhoton))->Et();
  tmva_ptom = tmva_pt/mass;
  
  tmva_badpf_iso = ((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad)*50/tmva_pt;
  tmva_goodpf_iso = ((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac)*50/tmva_pt;
  tmva_tkisopf = (*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt;
  //tmva_badpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad+0.0)*50/tmva_pt),10.);
  //tmva_goodpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac+0.0)*50/tmva_pt),10.);
  //tmva_tkisopf = min(sqrt((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt),10.);
  //tmva_badpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad+4.8)*50/tmva_pt),10.);
  //tmva_goodpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac+2.8)*50/tmva_pt),10.);
  //tmva_tkisopf = min(sqrt((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt),10.);
  if(MPDEBUG)  std::cout<<"BDT 4"<<std::endl;
  if(MPDEBUG)  std::cout<<"BDT 4"<<std::endl;
  
  tmva_sieie = ll->pho_sieie[iPhoton];
  tmva_drtotk = min(double(ll->pho_drtotk_25_99[iPhoton]),1.);
  tmva_drtotk = ll->pho_drtotk_25_99[iPhoton];
  tmva_hoe = ll->pho_hoe[iPhoton];
  tmva_r9 = ll->pho_r9[iPhoton];
  tmva_eta = fabs(((TLorentzVector*)ll->pho_p4->At(iPhoton))->Eta());
  
  if(MPDEBUG)  std::cout<<"BDT end"<<std::endl;
  Float_t mva = tmvaReader2->EvaluateMVA("Gradient");
  //cout<<"BDT  "<<jentry<<"  "<<iPhoton<<"  "<<vtx<<"  "<<isomax<<"  "<<badind<<"  "<<rhofac<<"  "<<rhofacbad<<"  "<<tmva_pt<<"  "<<tmva_badpf_iso<<"  "<<tmva_goodpf_iso<<"  "<<tmva_sieie<<"  "<<tmva_drtotk
  //    <<"  "<<tmva_hoe<<"  "<<tmva_tkisopf<<"  "<<tmva_r9<<"  "<<tmva_eta<<"  "<<mva<<endl;
  return mva;
}


Float_t StatAnalysisExclusive::BDT_ptom2(Int_t jentry, Int_t iPhoton, Int_t vtx, float mass) {
  // here we have a second BDT using ptom so we can compare developments more easily
  // we put in the diphoton mass just to scale the single photon pt (pt/m) 
  if(MPDEBUG)  std::cout<<"BDT"<<std::endl;

  int n_r9_categories = 3;
  int n_eta_categories = 2;
  tmva_cat = ll->PhotonCategory(iPhoton,n_r9_categories,n_eta_categories);

  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;
  ll->BdtGetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 after"<<std::endl;
  /*
  if(MPDEBUG)  std::cout<<"BDT 1"<<std::endl;

  if(jentry>0) {
  if (b_pho_sieie->GetReadEntry() != jentry)
  b_pho_sieie->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 1"<<std::endl;
  if (b_pho_pfiso_mycharged04->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 2"<<std::endl;
  if (b_pho_pfiso_myphoton04->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton04->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 3"<<std::endl;
  if (b_pho_pfiso_mycharged03->GetReadEntry() != jentry)
  b_pho_pfiso_mycharged03->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 4"<<std::endl;
  if (b_pho_pfiso_myphoton03->GetReadEntry() != jentry)
  b_pho_pfiso_myphoton03->GetEntry(jentry);
  //std::cout<<"BDT 1 - 5"<<std::endl;
  //if (b_pho_drtotk_25_99->GetReadEntry() != jentry)
  //b_pho_drtotk_25_99->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 6"<<std::endl;
  if (b_pho_hoe->GetReadEntry() != jentry)
  b_pho_hoe->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 7"<<std::endl;
  if (b_pho_r9->GetReadEntry() != jentry)
  b_pho_r9->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 8"<<std::endl;
  if (b_pho_p4->GetReadEntry() != jentry)
  b_pho_p4->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 9"<<std::endl;
  if (b_rho->GetReadEntry() != jentry)
  b_rho->GetEntry(jentry);
  }
  

  if(MPDEBUG)  std::cout<<"BDT 2"<<std::endl;
  */

  float isomax=-99;   int badind=0;
  for(int iv=0; iv<ll->vtx_std_n; iv++) if((*ll->pho_pfiso_mycharged04)[iPhoton][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged04)[iPhoton][iv]; }
  if(MPDEBUG)  std::cout<<"BDT 3"<<std::endl;
  
  float rhofacpf[6]={0.075, 0.082, 0.143, 0.050, 0.091, 0.106};          //move
  float rhofacbadpf[6]={0.141, 0.149, 0.208, 0.135, 0.162, 0.165};
  float rhofac=rhofacpf[tmva_cat];
  float rhofacbad=rhofacbadpf[tmva_cat];

  tmva_pt = ((TLorentzVector*)ll->pho_p4->At(iPhoton))->Et();
  tmva_ptom = tmva_pt/mass;
  
  //tmva_badpf_iso = ((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad)*50/tmva_pt;
  //tmva_goodpf_iso = ((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac)*50/tmva_pt;
  //tmva_tkisopf = (*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt;
  //tmva_badpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad+0.0)*50/tmva_pt),10.);
  //tmva_goodpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac+0.0)*50/tmva_pt),10.);
  //tmva_tkisopf = min(sqrt((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt),10.);
  tmva_badpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged04)[iPhoton][badind]+ll->pho_pfiso_myphoton04[iPhoton]-ll->rho*rhofacbad+4.8)*50/tmva_pt),10.);
  tmva_goodpf_iso = min(sqrt(((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]+ll->pho_pfiso_myphoton03[iPhoton]-ll->rho*rhofac+2.8)*50/tmva_pt),10.);
  tmva_tkisopf = min(sqrt((*ll->pho_pfiso_mycharged03)[iPhoton][vtx]*50/tmva_pt),(float) 10.);
  if(MPDEBUG)  std::cout<<"BDT 4"<<std::endl;
  if(MPDEBUG)  std::cout<<"BDT 4"<<std::endl;
  
  tmva_sieie = ll->pho_sieie[iPhoton];
  tmva_drtotk = min(double(ll->pho_drtotk_25_99[iPhoton]),1.);
  tmva_hoe = ll->pho_hoe[iPhoton];
  tmva_r9 = ll->pho_r9[iPhoton];
  tmva_eta = fabs(((TLorentzVector*)ll->pho_p4->At(iPhoton))->Eta());
  
  if(MPDEBUG)  std::cout<<"BDT end"<<std::endl;
  
  Float_t mva = tmvaReader3->EvaluateMVA("Gradient");
  //cout<<"BDT  "<<jentry<<"  "<<iPhoton<<"  "<<vtx<<"  "<<isomax<<"  "<<badind<<"  "<<rhofac<<"  "<<rhofacbad<<"  "<<tmva_pt<<"  "<<tmva_badpf_iso<<"  "<<tmva_goodpf_iso<<"  "<<tmva_sieie<<"  "<<tmva_drtotk
  //     <<"  "<<tmva_hoe<<"  "<<tmva_tkisopf<<"  "<<tmva_r9<<"  "<<tmva_eta<<"  "<<mva<<endl;
  return mva;
}


Float_t StatAnalysisExclusive::BDT_dipho(Int_t jentry, Int_t isl, Int_t il, float lmva, float slmva, float diphopt, float mass, float dmom) {


  ll->BdtGetEntry(jentry);

    /*
  if (b_pho_r9->GetReadEntry() != jentry)
  b_pho_r9->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 8"<<std::endl;
  if (b_pho_p4->GetReadEntry() != jentry)
  b_pho_p4->GetEntry(jentry);
    */
  if(MPDEBUG)  std::cout<<"BDT dipho"<<std::endl;

  tmva_leadr9 = ll->pho_r9[il];
  tmva_subleadr9 = ll->pho_r9[isl];
  tmva_leadeta = fabs(((TLorentzVector*)ll->pho_p4->At(il))->Eta());
  tmva_subleadeta = fabs(((TLorentzVector*)ll->pho_p4->At(isl))->Eta());
  tmva_etamax = max(fabs(((TLorentzVector*)ll->pho_p4->At(il))->Eta()),fabs(((TLorentzVector*)ll->pho_p4->At(isl))->Eta()));
  tmva_subleadptomass = ((TLorentzVector*)ll->pho_p4->At(isl))->Et()/mass;
  tmva_diphoptom = diphopt/mass;
  tmva_sumptom = (((TLorentzVector*)ll->pho_p4->At(il))->Et()+((TLorentzVector*)ll->pho_p4->At(isl))->Et())/mass;
  tmva_subleadmva = slmva;
  tmva_leadmva = lmva;
  tmva_dmom = dmom;
  
  Float_t mva = tmvaReader_dipho->EvaluateMVA("Gradient");
  
  return mva;
}

Float_t StatAnalysisExclusive::BDT_dipho2(Int_t jentry, Int_t isl, Int_t il, float lmva, float slmva, float diphopt, float mass) {

  ll->BdtGetEntry(jentry);

    /*
  if (b_pho_r9->GetReadEntry() != jentry)
  b_pho_r9->GetEntry(jentry);
  if(MPDEBUG)  std::cout<<"BDT 1 - 8"<<std::endl;
  if (b_pho_p4->GetReadEntry() != jentry)
  b_pho_p4->GetEntry(jentry);
    */

  if(MPDEBUG)  std::cout<<"BDT dipho"<<std::endl;

  tmva_leadr9 = ll->pho_r9[il];
  tmva_subleadr9 = ll->pho_r9[isl];
  tmva_leadeta = fabs(((TLorentzVector*)ll->pho_p4->At(il))->Eta());
  tmva_subleadeta = fabs(((TLorentzVector*)ll->pho_p4->At(isl))->Eta());
  tmva_etamax = max(fabs(((TLorentzVector*)ll->pho_p4->At(il))->Eta()),fabs(((TLorentzVector*)ll->pho_p4->At(isl))->Eta()));
  tmva_subleadptomass = ((TLorentzVector*)ll->pho_p4->At(isl))->Et()/mass;
  tmva_diphoptom = diphopt/mass;
  tmva_sumptom = (((TLorentzVector*)ll->pho_p4->At(il))->Et()+((TLorentzVector*)ll->pho_p4->At(isl))->Et())/mass;
  tmva_subleadmva = slmva;
  tmva_leadmva = lmva;
  
  Float_t mva = tmvaReader_dipho2->EvaluateMVA("Gradient");
  
  return mva;
}


void StatAnalysisExclusive::SetOutputNtupleVariables(int jentry, int itype, int leadind, int subleadind, int vtxind, float mass, TLorentzVector *leadp4, TLorentzVector *subleadp4, float evweight, float pileupWeight, TLorentzVector * jet1, TLorentzVector * jet2, TLorentzVector * jet3) {

  //FIX itype
  itype = itype_jim(itype);

  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables BEGIN"<<endl;
  //int leadind = dipho_leadind[diphoton_index];
  //int subleadind = dipho_subleadind[diphoton_index];

  //NOT NEEDED ANYMORE TLorentzVector newleadp4; //MARCO FIXED =  TLorentzVector(VertexCorrectedP4Hgg(leadind,vtxind));
  //NOT NEEDED ANYMORE TLorentzVector newsubleadp4;//MARCO FIXED  =  TLorentzVector(VertexCorrectedP4Hgg(subleadind,vtxind));

  //TVector3 leadcalopos = *((TVector3*)ll->pho_calopos->At(leadind));
  //TVector3 subleadcalopos = *((TVector3*)ll->pho_calopos->At(subleadind));



  t_j1pt=myAllLeadJPt;
  t_j1eta=myAllLeadJEta;
  t_j1phi=-20;
  if(jet1) t_j1phi=jet1->Phi();

  t_j2pt=myAllSubJPt;
  t_j2eta=myAllSubJEta;
  t_j2phi=-20;
  if(jet2) t_j2phi=jet2->Phi();

  t_jjmass=myAll_Mjj;
  t_jjdeta=myAlldEta;
  t_jjzep=myAllZep;
  t_jjdphi=myAlldPhi;


  t_j3pt=0.;
  t_j3eta=0.;
  t_j3phi=0.;

  //if(jet3)
	  //cout<<"AAA MARCOMM in set "<<jet3->Pt()<<endl;

  if(jet3) {
    t_j3pt=jet3->Pt();
    t_j3eta=jet3->Eta();
    t_j3phi=jet3->Phi();
  }

  Float_t leadetanonabs = (((TLorentzVector*)ll->sc_p4->At(ll->pho_scind[leadind]))->Eta());
  Float_t subleadetanonabs = (((TLorentzVector*)ll->sc_p4->At(ll->pho_scind[subleadind]))->Eta());
  Float_t leadeta = fabs(leadetanonabs);
  Float_t subleadeta = fabs(subleadetanonabs);

  TLorentzVector diphotonp4 = (*leadp4) + (*subleadp4);

  // apply k-factor to gammJet sample with 2 real photons. 
  // The k-factor from inputfiles is already in the weight, so here use kfactor = 1 except for the gammaJet 2-real case
  float k_factor = 1.;//mp->kfactor[indexfiles];
  //  if(itype==-1 || itype==-11 || itype==-21) {
  //    k_factor = GetWeightKfactor1D(HggKfactor1D, diphotonp4.M());
  //  } else if(itype==-2 || itype==-12 || itype==-22) {
  //    k_factor = GetWeightKfactor1D(HggKfactor1D, diphotonp4.M());
  //  } else if(itype==-3 && ll->pho_genmatch[leadind] && ll->pho_genmatch[subleadind]) {
  ////    k_factor = 1.7;
  //  }
  //std::cout << "itype: " << itype << "\tindexfiles: " << indexfiles << "\tkfac: " << k_factor << std::endl;

  t_run = ll->run;
  t_lumis = ll->lumis;
  t_event = ll->event;
  t_itype = itype;
  t_processid = ll->process_id;
  //   minuit_weight is the weight computed from inputfiles including the k-factor
  //   weight is not used to be myweight from this program which gives the PT-higgs correction and the efficiency correction for signal
  //cout<<"   weigth   "<<itype<<"  "<<minuit_weight<<"  "<<k_factor<<"  "<<pileup_reweight<<"  "<<t_w<<endl;

  //MARCO FIX higgs_genpt higgs_genpt higgs_genpt
  // if(jentry%1000==1&&itype>0&&itype%1000==1) cout<<jentry<<"  "<<itype<<"  mass="<<mass<<" out  higgs_genpt="<<ll->higgs_genpt<<"  diphopt="<<diphotonp4.Pt()<<"  t_w="<<t_w<<endl;
if(jentry%1000==1&&itype>0&&itype%1000==1) cout<<jentry<<"  "<<itype<<"  mass="<<mass<<" out  higgs_genpt= NOT THERE  diphopt="<<diphotonp4.Pt()<<"  t_w="<<t_w<<endl;



  //MARCO FIXED WEIGHT FIXED 
  //t_w = minuit_weight*k_factor*pileup_reweight*t_w;
  //t_wpu = pileup_reweight*weight;
  t_w =evweight;
  t_wpu = pileupWeight;


  t_mass = mass;
  t_dmom0=0.5*pow(pow(ll->pho_regr_energyerr[leadind]/ll->pho_regr_energy[leadind],2)+pow(ll->pho_regr_energyerr[subleadind]/ll->pho_regr_energy[subleadind],2),0.5);

  float leta=fabs( ((TLorentzVector*)ll->pho_p4->At(leadind))->Eta() );
  float seta=fabs( ((TLorentzVector*)ll->pho_p4->At(subleadind))->Eta() );

  //MARCO FIXED
  float leadErr = GetSmearSigma(leta,ll->pho_r9[leadind]);
  float subleadErr = GetSmearSigma(seta,ll->pho_r9[subleadind]);



  double errfrac=0.5*pow(pow(ll->pho_regr_energyerr[leadind]/ll->pho_regr_energy[leadind],2)+pow(leadErr,2)+pow(ll->pho_regr_energyerr[subleadind]/ll->pho_regr_energy[subleadind],2)+pow(subleadErr,2),0.5);
  t_dmom = errfrac;

  if(jentry%1000==1) {
    cout<<"jentry="<<jentry<<"  itype="<<itype;
    cout<<"  lres="<<ll->pho_regr_energyerr[leadind]/ll->pho_regr_energy[leadind]<<"  sres="<<ll->pho_regr_energyerr[subleadind]/ll->pho_regr_energy[subleadind]
	<<"  leta="<<leta<<"  lr9="<<ll->pho_r9[   leadind]<<"  "<<ll->PhotonCategory(leadind,2,2)<<"  lsmear="<<leadErr<<"  "
	<<"  seta="<<seta<<"  sr9="<<ll->pho_r9[subleadind]<<"  "<<ll->PhotonCategory(subleadind,2,2)<<"  ssmear="<<subleadErr;
    cout<<"    fractional error="<<errfrac<<endl;
  }

  t_deltaM=0;
  if(itype>0) t_deltaM=mass-float(itype/1000);

  t_diphocat2r92eta = ll->DiphotonCategory(leadind,subleadind,0.,2,2,1);

  //MARCO FIX CHECK
  //MARCO FIX NOW
  t_diphosubcat4 = diphoSubCategory(3,ll->DiphotonCategory(leadind,subleadind,0.,2,2,1),diphoCutLevel(leadind,subleadind,vtxind))-1;     //   jgb change replace chosen_vtx with vtxind (probably the same but its a parameter t this routine)

  t_category = ll->DiphotonCategory(leadind,subleadind,0.,4,3,1);
  t_barrel = (Int_t)(leadeta < 1.479 && subleadeta < 1.479);
  t_diphor9 = TMath::Min(ll->pho_r9[leadind],ll->pho_r9[subleadind]);
  t_diphoeta = fabs(diphotonp4.Eta());
  t_costhetastar = fabs(leadp4->P() - subleadp4->P())/diphotonp4.P();
  t_diphopt = diphotonp4.Pt();
  t_diphopz = diphotonp4.Pz();
  t_deltar = leadp4->DeltaR(*subleadp4);
  t_etamax = TMath::Max(leadeta,subleadeta);
  t_etamin = TMath::Min(leadeta,subleadeta);

  //MARCO CHANGED
  //t_deta = fabs(leadcalopos.Eta() - subleadcalopos.Eta());
  t_deta = fabs(leadetanonabs-subleadetanonabs);

  t_leadcat = ll->PhotonCategory(leadind,3,4);
  t_leadr9 = ll->pho_r9[leadind];
  //MARCO CHANGED
  //t_leadeta = leadcalopos.Eta();
  t_leadeta = leadetanonabs;


  t_leadpt = leadp4->Pt();
  t_leadgenmatch = GenIndexHgg(leadp4)>=0; //MARCO FIXED  = GenIndexHgg(leadp4)>=0;     
  t_leadbarrel = leadeta < 1.479;//ll->pho_barrel[leadind];
  t_leadhovere = ll->pho_hoe[leadind];
  t_leadsee = ll->pho_see[leadind];
  //t_leadspp = ll->pho_spp[leadind];
  t_leadpi0nn=0;

  t_leadecalhitsJ_060330; //MARCO FIX IGNORE = ll->pho_ecalJ_060330[leadind];

  t_leadtrkiso; //MARCO FIX IGNORE = ll->pho_ntrk_15_03[leadind];
  t_leadecaliso; //MARCO FIX IGNORE = ll->pho_ecalJ_060330[leadind];
  t_leadhcaliso; //MARCO FIX IGNORE = ll->pho_hcal_030[leadind];
  //t_leadtrkplusecal=t_leadtrkecone30 + t_leadecalhitsJ_060330;

  t_subleadcat = ll->PhotonCategory(subleadind,3,4);
  t_subleadr9 = ll->pho_r9[subleadind];
  //MARCO CHANGED
  //t_subleadeta = subleadcalopos.Eta();
  t_subleadeta = subleadetanonabs;

  t_subleadpt = subleadp4->Pt();
  //t_subleadgenmatch = ll->pho_genmatch[subleadind];
  t_subleadgenmatch = GenIndexHgg(subleadp4)>=0; //MARCO FIXED  = GenIndexHgg(subleadp4)>=0;     
  t_subleadbarrel = subleadeta < 1.479;//ll->pho_barrel[subleadind];
  t_subleadhovere = ll->pho_hoe[subleadind];
  t_subleadsee = ll->pho_see[subleadind];
  //t_subleadspp = ll->pho_spp[subleadind];
  t_subleadpi0nn=0;
  t_subleadecalhitsJ_060330; //MARCO FIX IGNORE = ll->pho_ecalJ_060330[subleadind];

  t_rho = ll->rho;
  t_nvtx = ll->vtx_std_n;



  if(!ll->GetCutValue("optree")) return;




  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 02  "<<endl;

  // cic index begin
//MARCO FIX CHECK
  t_leadci6cindex = PhotonCiCSelectionLevelJim(leadind,6,vtxind,0);
//MARCO FIX CHECK
  t_subleadci6cindex = PhotonCiCSelectionLevelJim(subleadind,6,vtxind,1);


  //MARCO FIX CHECK
  t_leadci6cpfindex = PhotonCiCpfSelectionLevelJim(leadind,6,vtxind,0);
  //MARCO FIX CHECK
  t_subleadci6cpfindex = PhotonCiCpfSelectionLevelJim(subleadind,6,vtxind,1);
  

  t_subleadci6cpfmva = 0.;
  t_subleadci6cpfmvacat = 0.;
  t_subleadci6cpfmvaptom = 0.;
  t_subleadci6cpfmvaptom2 = 0.;
  t_leadci6cpfmva = 0.;
  t_leadci6cpfmvacat = 0.;
  t_leadci6cpfmvaptom = 0.;
  t_leadci6cpfmvaptom2 = 0.;
  t_diphomva = 0.;
  t_diphomva2 = 0.;

  if(ll->GetCutValue("bdt")) {
    // FIXME 
    t_subleadci6cpfmva = 99;//BDT(jentry, subleadind,vtxind);
    t_subleadci6cpfmvacat = 99;//BDT_categorized(jentry, subleadind, vtxind, -1.);
    t_subleadci6cpfmvaptom = BDT_ptom(jentry, subleadind,vtxind, mass);
    t_subleadci6cpfmvaptom2 = 99;//BDT_ptom2(jentry, subleadind,vtxind, mass);
    t_leadci6cpfmva = 99;//BDT(jentry, leadind,vtxind);
    t_leadci6cpfmvacat = 99;//BDT_categorized(jentry, leadind,vtxind, 1.);
    t_leadci6cpfmvaptom = BDT_ptom(jentry, leadind,vtxind, mass);
    t_leadci6cpfmvaptom2 = 99;//BDT_ptom2(jentry, leadind,vtxind, mass);
    t_diphomva = BDT_dipho(jentry, leadind, subleadind, t_leadci6cpfmvaptom, t_subleadci6cpfmvaptom, t_diphopt, mass, t_dmom);
    t_diphomva2 = 99;//BDT_dipho2(jentry, leadind, subleadind, t_leadci6cpfmvaptom, t_subleadci6cpfmvaptom, t_diphopt, mass);
  }
  //std::pair<int,int> diphoton_indices(DiphotonCiCSelectionIndices( LEADCUTLEVEL, SUBLEADCUTLEVEL, leadPtMin, subleadPtMin, CICNCAT, -1, false));
  //leadind = diphoton_indices.first;
  //subleadind = diphoton_indices.second;

  t_leadcutindex = t_leadci6cpfindex;//PhotonCiCSelectionLevel(leadind,6,chosen_vtx,0);
  t_subleadcutindex = t_subleadci6cpfindex;//PhotonCiCSelectionLevel(subleadind,6,chosen_vtx,1);
//MARCO FIX CHECK
  t_leadci4cindex = PhotonCiCSelectionLevelJim(leadind,4,vtxind,0);
//MARCO FIX CHECK
  t_subleadci4cindex = PhotonCiCSelectionLevelJim(subleadind,4,vtxind,1);

  //std::pair<int,int> diphoton_indices(DiphotonCiCSelectionIndices( LEADCUTLEVEL, SUBLEADCUTLEVEL, leadPtMin, subleadPtMin, CICNCAT, -1, false));
  //leadind = diphoton_indices.first;
  //subleadind = diphoton_indices.second;
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 03  "<<endl;
  

  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 05  "<<endl;

  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 06  "<<endl;

  t_subleadtrkiso; //MARCO FIX IGNORE = ll->pho_ntrk_15_03[subleadind];
  t_subleadecaliso; //MARCO FIX IGNORE = ll->pho_ecalJ_060330[subleadind];
  t_subleadhcaliso; //MARCO FIX IGNORE = ll->pho_hcal_030[subleadind];

  //MARCO ????
  t_subleadtrkplusecal=t_subleadtrkecone30 + t_subleadecalhitsJ_060330;

  //new ntuple branches
  //  //lead vars
  t_leadpixel = ll->pho_haspixseed[leadind];
  t_leadsieie = ll->pho_sieie[leadind];
  t_leadtrkhollowdr03 = ll->pho_trksumpthollowconedr03[leadind];
  t_leadtrkhollowdr04 = ll->pho_trksumpthollowconedr04[leadind];
  t_leadtrksoliddr03 = ll->pho_trksumptsolidconedr03[leadind];
  t_leadtrksoliddr04 = ll->pho_trksumptsolidconedr04[leadind];
  t_leadecaldr03 = ll->pho_ecalsumetconedr03[leadind];
  t_leadecaldr04 = ll->pho_ecalsumetconedr04[leadind];
  t_leadhcaldr03 = ll->pho_hcalsumetconedr03[leadind];
  t_leadhcaldr04 = ll->pho_hcalsumetconedr04[leadind];
  //sublead vars
  t_subleadpixel = ll->pho_haspixseed[subleadind];
  t_subleadsieie = ll->pho_sieie[subleadind];
  t_subleadtrkhollowdr03 = ll->pho_trksumpthollowconedr03[subleadind];
  t_subleadtrkhollowdr04 = ll->pho_trksumpthollowconedr04[subleadind];
  t_subleadtrksoliddr03 = ll->pho_trksumptsolidconedr03[subleadind];
  t_subleadtrksoliddr04 = ll->pho_trksumptsolidconedr04[subleadind];
  t_subleadecaldr03 = ll->pho_ecalsumetconedr03[subleadind];
  t_subleadecaldr04 = ll->pho_ecalsumetconedr04[subleadind];
  t_subleadhcaldr03 = ll->pho_hcalsumetconedr03[subleadind];
  t_subleadhcaldr04 = ll->pho_hcalsumetconedr04[subleadind];
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 07  "<<endl;


  /*
  t_lead_tkiso_recvtx_030_005_0000_10_01    = ll->pho_tkiso_recvtx_030_005_0000_10_01[leadind];
  t_sublead_tkiso_recvtx_030_005_0000_10_01 = ll->pho_tkiso_recvtx_030_005_0000_10_01[subleadind];
  t_lead_tkiso_recvtx_030_006_0000_10_01    = ll->pho_tkiso_recvtx_030_006_0000_10_01[leadind];
  t_sublead_tkiso_recvtx_030_006_0000_10_01 = ll->pho_tkiso_recvtx_030_006_0000_10_01[subleadind];
  */

  t_lead_tkiso_recvtx_030_005_0000_10_01    = 100000.;
  t_sublead_tkiso_recvtx_030_005_0000_10_01 = 100000.;
  t_lead_tkiso_recvtx_030_006_0000_10_01    = 100000.;
  t_sublead_tkiso_recvtx_030_006_0000_10_01 = 100000.;
  

  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 14  "<<endl;




  // dxy
  /*
  t_lead_tkiso_recvtx_030_004_0000_10_02 = ll->pho_tkiso_recvtx_030_004_0000_10_02[leadind];
  t_sublead_tkiso_recvtx_030_004_0000_10_02 = ll->pho_tkiso_recvtx_030_004_0000_10_02[subleadind];
  */


  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 20  "<<endl;

  // vertex z
  //  gv_z = ((TVector3*)gv_pos->At(0))->Z();
  t_genvtxz; //MARCO FIX  LATER = ll->gv_z;
  t_recvtxz = ((TVector3*)ll->vtx_std_xyz->At(vtxind))->z();
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 21  "<<endl;

  //pfiso
  //cout<<"filling pfiso  "<<leadind<<"  "<<subleadind<<"  "<<vtxind<<"  "<<(void*)ll->pho_pfiso_mycharged03<<endl;
  //cout<<(*ll->pho_pfiso_mycharged03)[0][0]<<"  "<<ll->pho_pfiso_mycharged03->size()<<" "<<(*ll->pho_pfiso_mycharged03)[0].size()<<endl;

  t_lead_pfiso_charged03 = (*ll->pho_pfiso_mycharged03)[leadind][vtxind];                               //   jgb    _noveto  04
  t_lead_pfiso_photon03 = ll->pho_pfiso_myphoton03[leadind];                               //   jgb    _noveto  04
  t_lead_pfiso_neutral03 = ll->pho_pfiso_myneutral03[leadind];                               //   jgb    _noveto  04
  t_lead_pfiso_charged04 = (*ll->pho_pfiso_mycharged04)[leadind][vtxind];                               //   jgb    _noveto  04
  t_lead_pfiso_photon04 = ll->pho_pfiso_myphoton04[leadind];                               //   jgb    _noveto  04
  t_lead_pfiso_neutral04 = ll->pho_pfiso_myneutral04[leadind];                               //   jgb    _noveto  04
  t_sublead_pfiso_charged03 = (*ll->pho_pfiso_mycharged03)[subleadind][vtxind];                               //   jgb    _noveto  04
  t_sublead_pfiso_photon03 = ll->pho_pfiso_myphoton03[subleadind];                               //   jgb    _noveto  04
  t_sublead_pfiso_neutral03 = ll->pho_pfiso_myneutral03[subleadind];                               //   jgb    _noveto  04
  t_sublead_pfiso_charged04 = (*ll->pho_pfiso_mycharged04)[subleadind][vtxind];                               //   jgb    _noveto  04
  t_sublead_pfiso_photon04 = ll->pho_pfiso_myphoton04[subleadind];                               //   jgb    _noveto  04
  t_sublead_pfiso_neutral04 = ll->pho_pfiso_myneutral04[subleadind];                               //   jgb    _noveto  04
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 22  "<<endl;

  float isomax=-99, dzmin=99;  int badind=0, genind=0;
  for(int iv=0; iv<t_nvtx; iv++) if(fabs(t_genvtxz-((TVector3*)ll->vtx_std_xyz->At(iv))->z())<dzmin) {genind=iv; dzmin=fabs(t_genvtxz-((TVector3*)ll->vtx_std_xyz->At(iv))->z()); }
  for(int iv=0; iv<t_nvtx; iv++) if((*ll->pho_pfiso_mycharged03)[leadind][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged03)[leadind][iv]; }
  t_lead_pfiso_charged_badvtx_03 = (*ll->pho_pfiso_mycharged03)[leadind][badind];                               //   jgb    _noveto  04
  t_lead_pfiso_photon_badvtx_03 = ll->pho_pfiso_myphoton03[leadind];                               //   jgb    _noveto  04
  t_lead_pfiso_neutral_badvtx_03 = ll->pho_pfiso_myneutral03[leadind];                               //   jgb    _noveto  04
  isomax=-99;  badind=0;
  for(int iv=0; iv<t_nvtx; iv++) if((*ll->pho_pfiso_mycharged04)[leadind][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged04)[leadind][iv]; }
  t_lead_pfiso_charged_badvtx_04 = (*ll->pho_pfiso_mycharged04)[leadind][badind];                               //   jgb    _noveto  04
  t_lead_pfiso_photon_badvtx_04 = ll->pho_pfiso_myphoton04[leadind];                               //   jgb    _noveto  04
  t_lead_pfiso_neutral_badvtx_04 = ll->pho_pfiso_myneutral04[leadind];                               //   jgb    _noveto  04
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 23  "<<endl;

  isomax=-99;  badind=0;
  for(int iv=0; iv<t_nvtx; iv++) if((*ll->pho_pfiso_mycharged04)[subleadind][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged04)[subleadind][iv]; }
  t_sublead_pfiso_charged_badvtx_04 = (*ll->pho_pfiso_mycharged04)[subleadind][badind];                               //   jgb    _noveto  04
  t_sublead_pfiso_photon_badvtx_04 = ll->pho_pfiso_myphoton04[subleadind];                               //   jgb    _noveto  04
  t_sublead_pfiso_neutral_badvtx_04 = ll->pho_pfiso_myneutral04[subleadind];                               //   jgb    _noveto  04
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 24  "<<endl;

  isomax=-99;  badind=0;
  for(int iv=0; iv<t_nvtx; iv++) if((*ll->pho_pfiso_mycharged03)[subleadind][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged03)[subleadind][iv]; }
  t_sublead_pfiso_charged_badvtx_03 = (*ll->pho_pfiso_mycharged03)[subleadind][badind];                               //   jgb    _noveto  04
  t_sublead_pfiso_photon_badvtx_03 = ll->pho_pfiso_myphoton03[subleadind];                               //   jgb    _noveto  04
  t_sublead_pfiso_neutral_badvtx_03 = ll->pho_pfiso_myneutral03[subleadind];                               //   jgb    _noveto  04
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 25  "<<endl;

  t_badvtxz = ((TVector3*)ll->vtx_std_xyz->At(badind))->z();
  /*
  if((itype<-9&&!t_subleadgenmatch&&fabs(t_recvtxz-t_genvtxz)>1||itype>0&&t_subleadgenmatch&&fabs(t_recvtxz-t_genvtxz)<1)&&t_sublead_pfiso_charged_badvtx_03-t_sublead_pfiso_charged03>3&&t_sublead_pfiso_charged_badvtx_03>9) {
    cout<<"itype="<<itype<<"  diphopt="<<t_diphopt<<"  t_genvtxz="<<t_genvtxz<<"  recvtx="<<vtxind<<"  "<<t_recvtxz<<"  badvtx="<<badind<<"  "<<t_badvtxz
	<<"  pfisoch="<<t_sublead_pfiso_charged03<<"  pfisoch_bad="<<t_sublead_pfiso_charged_badvtx_03<<endl;
    for(int v=0; v<vtx_std_n; v++) cout<<vtx_std_scalarpt[v]<<"  "; cout<<endl;
  }
  */
  t_genvtx_sumpt=ll->vtx_std_scalarpt[genind];
  t_recvtx_sumpt=ll->vtx_std_scalarpt[vtxind];
  t_badvtx_sumpt=ll->vtx_std_scalarpt[badind];
  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables 26  "<<endl;

  if(jentry%1000==1) {
    cout<<subleadind<<" "<<leadind
	<<"  vtxind="<<vtxind<<"  badind="<<badind
	<<"  t_subleadsieie "<<t_subleadsieie
	<<"  t_subleadhovere "<<t_subleadhovere
	<<"  t_subleadr9 "<<t_subleadr9
	<<"  t_sublead_pfiso_charged03 "<<t_sublead_pfiso_charged03
	<<"  t_sublead_pfiso_photon03 "<<t_sublead_pfiso_photon03
	<<"  t_sublead_pfiso_charged_badvtx_03 "<<t_sublead_pfiso_charged_badvtx_03
	<<"   t_sublead_pfiso_neutral "<<t_sublead_pfiso_neutral04<<endl;
    //	<<"t_sublead_drtotk_25_99 "<<t_subleaddrtotk_25_99
  }


  t_lead_drtotk_25_99 = ll->pho_drtotk_25_99[leadind];
  t_sublead_drtotk_25_99 = ll->pho_drtotk_25_99[subleadind];

  //  pf first guess isolation variables   pfiso


  if((optree->GetEntries()<100||optree->GetEntries()%1000==1)&&t_mass>95&&t_mass<145&&t_leadci6cindex>2&&t_subleadci6cindex>2) {
    float stkreciso=(*(ll->pho_tkiso_recvtx_030_002_0000_10_01))[subleadind][vtxind];
    float ltkreciso=(*(ll->pho_tkiso_recvtx_030_002_0000_10_01))[leadind][vtxind];

//MARCO FIX CHECK

    cout<<"setting optree variables run "
	<<ll->run<<" ev  "<<ll->event<<" cat "<<t_diphocat2r92eta<<" dscat "<<t_diphosubcat4<<" dcutlev "<<diphoCutLevel(leadind,subleadind,vtxind)<<" lind,sind,vtxind "<<leadind<<subleadind<<vtxind<<" lcut  "<<t_leadci6cindex<<" scut "<<t_subleadci6cindex<<" ltkreciso "
	<<ltkreciso<<" stkrecios "<<stkreciso<<" mass "<<t_mass<<" w "<<t_w<<"  diphomva="<<t_diphomva<<endl;
  }

  if(MPDEBUGCOPY) cout << "SetOutputNtupleVariables END"<<endl;
}

std::vector<Float_t> StatAnalysisExclusive::arrayToVector(size_t length, const Float_t *data) {
  std::vector<float> retval;
  std::copy (data, data + length,std::back_inserter ( retval ) );
  return retval;
}
float StatAnalysisExclusive::GetSmearSigma(float eta, float r9, int epoch){
  // not using epoch for now
  // r9cat + nr9Cat*etaCat
  float sigma = -1;
  if(epoch==0) {
    int nr9Cat=2;
    int r9Cat=(int)(r9<0.94);
    int nEtaCat=4;
    int etaCat=(int)(eta>1.0) + (int)(eta>1.5) + (int)(eta>2.0);
    sigma=smear_nov14[(int)(r9Cat + nr9Cat*etaCat)];
  }

  return sigma;
}



Int_t StatAnalysisExclusive::GenIndexHgg(TLorentzVector *p_p4) {

  /*
    gp_n
    gp_p4
    gp_status
    gp_pdgid
    gp_mother

   Int_t           gp_n;
   TClonesArray    *gp_p4;
   ////TClonesArray    *gp_vtx;
   Short_t         gp_status[1538];   //[gp_n]
   Short_t         gp_pdgid[1538];   //[gp_n]
   Short_t         gp_mother[1538];   //[gp_n]
  */

  Float_t deltaRmin = 100;
  Int_t genindex = -1;
  //std::cout << "gp_n: " << ll->gp_n << std::endl;
  for(Int_t i=0; i!=ll->gp_n; ++i) {
    //std::cout << "i: " << i << std::endl;
    if(ll->gp_status[i]!=1) continue;
    if(fabs(ll->gp_pdgid[i])!=22) continue; 
    if(ll->gp_status[ll->gp_mother[i]]!=3) continue;
    TLorentzVector * genp4 = (TLorentzVector *) ll->gp_p4->At(i);
    //std::cout << "ll->gp_pdgid[i]: " << ll->gp_pdgid[i] << "gp_status[i]: " << ll->gp_status[i] << "gp_mother[i]: " << ll->gp_mother[i] << "gp_status["<<ll->gp_mother[i]<<"]: " << ll->gp_status[ll->gp_mother[i]] << "gp_pdgid["<<ll->gp_mother[i]<<"]: " << ll->gp_pdgid[ll->gp_mother[i]] << "genp4->Et: " << genp4->Et() << std::endl;
    //std::cout << "gp_status[i]: " << ll->gp_status[i] << std::endl;
    //std::cout << "gp_mother[i]: " << ll->gp_mother[i] << std::endl;
    //std::cout << "gp_status["<<ll->gp_mother[i]<<"]: " << ll->gp_status[ll->gp_mother[i]] << std::endl;
    Float_t deltaR = genp4->DeltaR(*p_p4);
    //std::cout << "deltaR : " << deltaR  << std::endl;
    if( deltaR < deltaRmin ) {
      genindex=i;
      deltaRmin = deltaR;
      if(deltaRmin < 0.08) return genindex;
    }
  }
  if(deltaRmin > 0.06) return -1;
  return genindex;
}


Int_t StatAnalysisExclusive::itype_jim(int itype) {

  if(itype==0) {
    return 0;
  }
  else if(itype>0) {
    if(itype==1) return -24; //qcdff 30
    //if(itype==2) return -24; //qcdff 40
    if(itype==3) return -13; //gj_pf 20
    if(itype==4) return -1;  //diphojets
    if(itype==5) return -2;  //box25-250
    if(itype==6) return -25; //DY
    //if(itype==7) return ;
    if(itype==8) return -2;  //box25-250
    //if(itype==9) return ;
    //if(itype==10) return ;
    if(itype==11) return -14; //qcdpf 30
    if(itype==12) return -14; //qcdpf 40
    if(itype==2) return -14; //qcdpf 40

  }
  else if(itype<0) {

  //MH=100 -69

  if(itype==-69) return 100001;
  if(itype==-70) return 100002;
  if(itype==-72) return 100003;
  if(itype==-71) return 100004;

  //MH=105 -13

  if(itype==-13) return 105001;
  if(itype==-14) return 105002;
  if(itype==-16) return 105003;
  if(itype==-15) return 105004;

  //MH=110 -17

  if(itype==-17) return 110001;
  if(itype==-18) return 110002;
  if(itype==-20) return 110003;
  if(itype==-19) return 110004;

  //MH=115 -21

  if(itype==-21) return 115001;
  if(itype==-22) return 115002;
  if(itype==-24) return 115003;
  if(itype==-23) return 115004;

  //MH=120 -25

  if(itype==-25) return 120001;
  if(itype==-26) return 120002;
  if(itype==-28) return 120003;
  if(itype==-27) return 120004;

  //MH=121 -53

  if(itype==-53) return 121001;
  if(itype==-54) return 121002;
  if(itype==-56) return 121003;
  if(itype==-55) return 121004;

  //MH=123 -57

  if(itype==-57) return 123001;
  if(itype==-58) return 123002;
  if(itype==-60) return 123003;
  if(itype==-59) return 123004;

  //MH=125 -37

  if(itype==-37) return 125001;
  if(itype==-38) return 125002;
  if(itype==-40) return 125003;
  if(itype==-39) return 125004;

  //MH=130 -29

  if(itype==-29) return 130001;
  if(itype==-30) return 130002;
  if(itype==-32) return 130003;
  if(itype==-31) return 130004;
  //MH=135 41

  if(itype==-41) return 135001;
  if(itype==-42) return 135002;
  if(itype==-44) return 135003;
  if(itype==-43) return 135004;

  //MH=140 33

  if(itype==-33) return 140001;
  if(itype==-34) return 140002;
  if(itype==-36) return 140003;
  if(itype==-35) return 140004;

  //MH=145 45

  if(itype==-45) return 145001;
  if(itype==-46) return 145002;
  if(itype==-48) return 145003;
  if(itype==-47) return 145004;

  //MH=150 49
  if(itype==-49) return 150001;
  if(itype==-50) return 150002;
  if(itype==-52) return 150003;
  if(itype==-51) return 150004;
  }

  std::cout<<" ATTENTION WRONG ITYPE: "<<itype<<std::endl;
  return 99999;

}




// ---------------------------------------------------------------------------------------------------------------------------------------------
std::pair<int,int> StatAnalysisExclusive::DiphotonCiCSelectionIndicesJim( phoCiCIDLevel LEADCUTLEVEL, phoCiCIDLevel SUBLEADCUTLEVEL, Float_t leadPtMin, Float_t subleadPtMin, int ncat, int vtxind, bool applyPtoverM) {

  if(MPDEBUGCOPY) {
    cout<<"LEADCUTLEVEL / SUBLEADCUTLEVEL / leadPtMin / subleadPtMin / ncat / vtxind / applyPtoverM"<<endl;
    cout<<LEADCUTLEVEL<<"\t"<<SUBLEADCUTLEVEL<<"\t"<<leadPtMin<<"\t"<<subleadPtMin<<"\t"<<ncat<<"\t"<<vtxind<<"\t"<<applyPtoverM<<endl;
  }


  //mmmmmmmm
  int myprint = 0;

  int selected_lead_index = -1;
  int selected_sublead_index = -1;
  float selected_lead_pt = -1;
  float selected_sublead_pt = -1;

  float lim1=1.4442;
  float lim2=1.566;
  float lim3=2.5;

  int diphoind=-1;
  for(int ipho=0;ipho!=ll->pho_n;++ipho) {
    for(int iipho=0;iipho!=ll->pho_n;++iipho) {
      if(iipho == ipho)continue;


      if(MPDEBUGCOPY) 
	cout<<"looking1 for dipho ind: "<<vtxind<<endl;
      if(vtxind<0) {
	//get right index to photons and variables....
	for(int idipho=0;idipho!=ll->dipho_n;++idipho) {
	  if(MPDEBUGCOPY) cout<<"looking for dipho ind: ipho, dipho"<<ipho<<" "<<iipho<<" "<<ll->dipho_leadind[idipho]<<" "<<ll->dipho_subleadind[idipho]<<endl;
	  if(ipho==ll->dipho_leadind[idipho]&&iipho==ll->dipho_subleadind[idipho]) {
	    vtxind=ll->dipho_vtxind[idipho];                 //   bug??   we should not reset the input parameter vtxind in this loop.  Should work with some copy.
	    diphoind=idipho;
	  }
	}
      }
      else {
	for(int idipho=0;idipho!=ll->dipho_n;++idipho) {
	  if(MPDEBUGCOPY) cout<<"looking2 for dipho ind: ipho, dipho"<<ipho<<" "<<iipho<<" "<<ll->dipho_leadind[idipho]<<" "<<ll->dipho_subleadind[idipho]<<endl;
	  if(ipho==ll->dipho_leadind[idipho]&&iipho==ll->dipho_subleadind[idipho]) {
	    diphoind=idipho;
	  }
	}
	if(diphoind<0) cout<<" **** WARNING: No diphot ind found for the requested vtxind: "<<vtxind<<endl;

      }

      if(MPDEBUGCOPY) cout<<"here marco11"<<ipho<<" "<<vtxind<<endl;


      //TLorentzVector * iphop4 = (TLorentzVector*)pho_p4->At(ipho);
      //TLorentzVector iphop4; //MARCO FIX  = TLorentzVector(VertexCorrectedP4Hgg(ipho,vtxind));
      //MARCO FIX CHECK
      TLorentzVector iphop4 = ll->get_pho_p4(ipho,vtxind);


      if(MPDEBUGCOPY) cout<<"here marco11"<<endl;
      //mmmmmm problem sean
      //if(iphop4.Et() < leadPtMin || fabs(iphop4.Eta()) > 2.5)continue;
      if(iphop4.Et() < leadPtMin) continue;
      if(MPDEBUGCOPY) cout<<"here marco12"<<endl;
      if(iphop4.Et() < selected_lead_pt) continue;
      if(MPDEBUGCOPY) cout<<"here marco13"<<endl;
      TVector3 * caloposi = (TVector3*)ll->pho_calopos->At(ipho);
      if(MPDEBUGCOPY) cout<<"here marco14"<<endl;
      Float_t etai = fabs(((TLorentzVector*)ll->sc_p4->At(ll->pho_scind[ipho]))->Eta());
      if(MPDEBUGCOPY) cout<<"here marco15"<<endl;

      if(fabs(etai)>lim1&&fabs(etai)<lim2) continue;
      if(MPDEBUGCOPY) cout<<"here marco16"<<endl;
      if(fabs(etai)>lim3) continue;
      if(MPDEBUGCOPY) cout<<"here marco17"<<endl;
      if(MPDEBUGCOPY) cout<<"here marco18"<<endl;

      if(PhotonCiCSelectionLevelJim(ipho, ncat, vtxind, 0, diphoind, myprint) < LEADCUTLEVEL) continue;

      //TLorentzVector * iiphop4 = (TLorentzVector*)ll->pho_p4->At(iipho);
      //mmmmmm problem sean
      //TLorentzVector iiphop4; //MARCO FIX  = TLorentzVector(VertexCorrectedP4Hgg(iipho,vtxind));
      //MARCO FIX CHECK
      TLorentzVector iiphop4 = ll->get_pho_p4(iipho,vtxind);



      //if(iiphop4.Et() < subleadPtMin || fabs(iiphop4.Eta()) > 2.5)continue;
      if(iiphop4.Et() >= iphop4.Et())continue;
      if(iiphop4.Et() < subleadPtMin) continue;
      if(iiphop4.Et() < selected_sublead_pt) continue;
      TVector3 * caloposii = (TVector3*)ll->pho_calopos->At(iipho);
      Float_t etaii = fabs(((TLorentzVector*)ll->sc_p4->At(ll->pho_scind[iipho]))->Eta());

      if(MPDEBUGCOPY) cout<<"here marco13"<<endl;

      if(fabs(etaii)>lim1&&fabs(etaii)<lim2) continue;
      if(fabs(etaii)>lim3) continue;

      float m_gamgam = (iphop4+iiphop4).M();
      float L_ptom = iphop4.Et()/m_gamgam;
      float S_ptom = iiphop4.Et()/m_gamgam;
      if(applyPtoverM && (L_ptom < 0.33 || S_ptom<0.25)) continue;
      if(MPDEBUGCOPY) cout<<"here marco14"<<endl;

      if(PhotonCiCSelectionLevelJim(iipho, ncat, vtxind, 1, diphoind, myprint) < SUBLEADCUTLEVEL) continue;
      if(MPDEBUGCOPY) cout<<"here marco15"<<endl;
      // if here, diphoton passed all cuts.
      //std::cout << "FOUND DIPHOTON" << std::endl;
      //
      //if(iphop4.Et() > selected_lead_pt) //READDED MARCO 23/7/11
      {
        selected_lead_pt = iphop4.Et();
        selected_sublead_pt = iiphop4.Et();
        selected_lead_index = ipho;
        selected_sublead_index = iipho;
      }

    }// end photon loop (iipho), aka sublead
  }// end photon loop (ipho), aka lead
      if(MPDEBUGCOPY) cout<<"here marco16"<<endl;

  std::pair<int,int> dipho_inds(selected_lead_index,selected_sublead_index);
  return dipho_inds;
      if(MPDEBUGCOPY) cout<<"here marco17"<<endl;

}


// ---------------------------------------------------------------------------------------------------------------------------------------------
int   StatAnalysisExclusive::PhotonCiCSelectionLevelJim( int photon_index, int ncat, int vtxind, int doSublead, int diphoind, int print) {
  //  This routine picks out the correct values of the isolation etc. so that the photon id can be applied.
 
  int cutlevelpassed = -1;

  int n_r9_categories = -1;
  int n_eta_categories = -1;
  if(ncat==6) {
    n_r9_categories = 3;
    n_eta_categories = 2;
  } else if(ncat==4) {
    n_r9_categories = 2;
    n_eta_categories = 2;
  }
  int photon_category = ll->PhotonCategory(photon_index,n_r9_categories,n_eta_categories);

  //TLorentzVector * phop4 = (TLorentzVector*)ll->pho_p4->At(photon_index);
  //TLorentzVector phop4; //MARCO FIX  = TLorentzVector(VertexCorrectedP4Hgg(photon_index,vtxind));
  //MARCO FIX CHECK
  TLorentzVector phop4 = ll->get_pho_p4(photon_index,vtxind);



  // MARCO FIX float val_tkiso = ll->pho_tkiso_recvtx_030_002_0000_10_01[photon_index];
  float val_tkiso = (*(ll->pho_tkiso_recvtx_030_002_0000_10_01))[photon_index][vtxind];
  float val_ecaliso = ll->pho_ecalsumetconedr03[photon_index];
  float val_hcaliso = ll->pho_hcalsumetconedr04[photon_index];
  float val_ecalisobad = ll->pho_ecalsumetconedr04[photon_index];
  float val_hcalisobad = ll->pho_hcalsumetconedr04[photon_index];
  float val_tkisobad = ll->pho_tkiso_badvtx_040_002_0000_10_01[photon_index];
  float val_drtotk_25_99 = ll->pho_drtotk_25_99[photon_index];
  float	val_pfiso_charged = (*ll->pho_pfiso_mycharged03)[photon_index][vtxind];                               //   jgb    _noveto  04
  float	val_pfiso_photon = ll->pho_pfiso_myphoton03[photon_index];                               //   jgb    _noveto  04
  float	val_pfiso_neutral = ll->pho_pfiso_myneutral03[photon_index];                               //   jgb    _noveto  04

  if(MPDEBUG) cout<<"DRTOTK: "<<photon_index<<" "<<diphoind<<" "<<ll->pho_drtotk_25_99[photon_index]<<" "
		  //MARCO FIX CHECK << ll->dipho_lead_pho_drtotk_25_99[diphoind]<<" "
		  //MARCO FIX CHECK <<ll->dipho_sublead_pho_drtotk_25_99[diphoind] <<" "
		  << ll->pho_drtotk_25_99[ll->dipho_leadind[diphoind]]<<" "
		  << ll->pho_drtotk_25_99[ll->dipho_subleadind[diphoind]]<<" "
		  <<endl;


  //????
  // MARCO FIX CHECK WELL I CHANGED IT
  //????
  //????

  /*
  if(vtxind>=0) {
    if(diphoind!=-1) {
      if(!doSublead) {
	 //MARCO FIX CHECK val_drtotk_25_99; //MARCO FIX !!!!! simple  = ll->dipho_lead_pho_drtotk_25_99[diphoind];
	//MARCO FIX CHECK val_tkiso; //MARCO FIX  = ll->dipho_lead_pho_tkiso_recvtx_030_002_0000_10_01[diphoind];
	//MARCO FIX CHECK val_tkisobad; //MARCO FIX = ll->dipho_lead_pho_tkiso_badvtx_040_002_0000_10_01[diphoind]; 
	val_drtotk_25_99 = ll->pho_drtotk_25_99[dipho_leadind[diphoind]];
	val_tkiso = ll->pho_tkiso_recvtx_030_002_0000_10_01[dipho_leadind[diphoind]];
	val_tkisobad = ll->pho_tkiso_badvtx_040_002_0000_10_01[dipho_leadind[diphoind]]; 
	val_pfiso_charged = (*ll->pho_pfiso_mycharged03)[photon_index][dipho_leadind[diphoind]];                               //   jgb    _noveto  04
	val_pfiso_photon = ll->pho_pfiso_myphoton03[photon_index];                               //   jgb    _noveto  04
	val_pfiso_neutral = ll->pho_pfiso_myneutral03[photon_index];                               //   jgb    _noveto  04

      }
      else {
	val_drtotk_25_99; //MARCO FIX = ll->dipho_sublead_pho_drtotk_25_99[diphoind];
	val_tkiso; //MARCO FIX = ll->dipho_sublead_pho_tkiso_recvtx_030_002_0000_10_01[diphoind];
	val_tkisobad; //MARCO FIX = ll->dipho_sublead_pho_tkiso_badvtx_040_002_0000_10_01[diphoind];
      }
    }
  }
*/
  val_drtotk_25_99 = ll->pho_drtotk_25_99[photon_index];
  //val_tkiso = ll->pho_tkiso_recvtx_030_002_0000_10_01[photon_index];
  val_tkiso = (*(ll->pho_tkiso_recvtx_030_002_0000_10_01))[photon_index][vtxind];

  val_tkisobad = ll->pho_tkiso_badvtx_040_002_0000_10_01[photon_index]; 
  val_pfiso_charged = (*ll->pho_pfiso_mycharged03)[photon_index][vtxind];                               //   jgb    _noveto  04
  val_pfiso_photon = ll->pho_pfiso_myphoton03[photon_index];                               //   jgb    _noveto  04
  val_pfiso_neutral = ll->pho_pfiso_myneutral03[photon_index];                               //   jgb    _noveto  04


  //????
  // END MARCO FIX CHECK WELL I CHANGED IT
  //????


  float val_sieie = ll->pho_sieie[photon_index];
  float val_hoe = ll->pho_hoe[photon_index];
  float val_r9 = ll->pho_r9[photon_index];
  float val_pixel = (float)ll->pho_haspixseed[photon_index];







  //float rhofacbad=0.40, rhofac=0.05;
  float rhofacbad=0.52, rhofac=0.17;
  float val_isosum=(val_tkiso+val_ecaliso+val_hcaliso-ll->rho*rhofac);
  float val_isosumbad=(val_tkisobad+val_ecalisobad+val_hcalisobad-ll->rho*rhofacbad);
  float val_trkiso=(val_tkiso);
  float val_et=phop4.Et();

  //mmmmmmmmmm
  //1 1 69613 7.6 GeV isolation

  if(print) {
    cout<<photon_index<<" "<<endl;
    cout<<"val_isosumoet "<<(val_tkiso+val_ecaliso+val_hcaliso-ll->rho*rhofac)<<endl;
    cout<<"val_isosumoetbad "<<(val_tkisobad+val_ecalisobad+val_hcalisobad-ll->rho*rhofacbad)<<endl;
    cout<<"val_trkisooet "<<(val_tkiso)<<endl;
    cout<<"val_sieie "<<ll->pho_sieie[photon_index]<<endl;
    cout<<"val_hoe "<<ll->pho_hoe[photon_index]<<endl;
    cout<<"val_r9 "<<ll->pho_r9[photon_index]<<endl;
    cout<<"val_drtotk_25_99 "<<ll->pho_drtotk_25_99[photon_index]<<endl;
    //cout<<"val_pixel "<<(float)ll->pho_haspixseed[photon_index]<<endl;
    
    
    //MARCO FIX cout<<"              val_tkiso "<<ll->pho_tkiso_recvtx_030_002_0000_10_01[photon_index]<<endl;
  cout<<"              val_tkiso "<<(*(ll->pho_tkiso_recvtx_030_002_0000_10_01))[photon_index][vtxind]<<endl;



    cout<<"              val_ecaliso "<<ll->pho_ecalsumetconedr03[photon_index]<<endl;
    cout<<"              val_hcaliso "<<ll->pho_hcalsumetconedr04[photon_index]<<endl;
    cout<<"              val_ecalisobad "<<ll->pho_ecalsumetconedr04[photon_index]<<endl;
    cout<<"              val_hcalisobad "<<ll->pho_hcalsumetconedr04[photon_index]<<endl;
    cout<<"              val_tkisobad "<<ll->pho_tkiso_badvtx_040_002_0000_10_01[photon_index]<<endl;
    
  }
  float val_eta = fabs(((TVector3*)ll->pho_calopos->At(photon_index))->Eta());

  //  Here we use the CiC routines to calculate the photon cut level for either 4 or 6 categories
  //  The parameters are the cut variables so that the choice of the photons, vertex etc. can be separated from the basic Id

  //cout<<"   PhotonCiCSelectionLevel  "<<photon_index<<"  "<<ncat<<"  "<<vtxind<<"  "<<doSublead<<"  "<<diphoind<<"  "<<val_r9<<"  "<<val_eta<<"  "<<val_isosumoet<<"  "
  //    <<val_isosumoetbad<<"  "<<val_trkisooet<<"  "<<val_sieie<<"  "<<val_hoe<<"  "<<val_drtotk_25_99<<endl;

  if(ncat==4) { cutlevelpassed = photonCutLevel4(val_r9, val_eta, val_et, val_isosum, val_isosumbad, val_trkiso, val_sieie, val_hoe, val_drtotk_25_99); }
  else if(ncat==6) { cutlevelpassed = photonCutLevel6(val_r9, val_eta, val_et, val_isosum, val_isosumbad, val_trkiso, val_sieie, val_hoe, val_drtotk_25_99);}
  else { cerr<<"Photon selection for "<<ncat<<" categories does not exist"<<endl; }


  return cutlevelpassed;

}



// ---------------------------------------------------------------------------------------------------------------------------------------------
int   StatAnalysisExclusive::PhotonCiCpfSelectionLevelJim( int photon_index, int ncat, int vtxind, int doSublead, int diphoind, int print) {
  //  This routine picks out the correct values of the isolation etc. so that the photon id can be applied.
  //  This version is for the particle flow isolation
  if(MPDEBUGCOPY) cout << "PhotonCiCpfSelectionLevel  BEGIN"<<endl;
 
  int cutlevelpassed = -1;

  int n_r9_categories = 3;
  int n_eta_categories = 2;
  int photon_category = ll->PhotonCategory(photon_index,n_r9_categories,n_eta_categories);

  //TLorentzVector phop4; //MARCO FIX  = TLorentzVector(VertexCorrectedP4Hgg(photon_index,vtxind));
  //MARCO FIX CHECK
  TLorentzVector phop4 = ll->get_pho_p4(photon_index,vtxind);


  float val_et=phop4.Et();
  float val_eta = fabs(((TVector3*)ll->pho_calopos->At(photon_index))->Eta());

  if(MPDEBUGCOPY) cout << "PhotonCiCpfSelectionLevel  01"<<endl;
  //pfiso
  //cout<<"filling pfiso  "<<leadind<<"  "<<subleadind<<"  "<<vtxind<<"  "<<(void*)>ll->pho_pfiso_mycharged03<<endl;
  //cout<<(*ll->pho_pfiso_mycharged03)[0][0]<<"  "<<ll->pho_pfiso_mycharged03->size()<<" "<<(*ll->pho_pfiso_mycharged03)[0].size()<<endl;
  float val_pfiso_charged03 = (*ll->pho_pfiso_mycharged03)[photon_index][vtxind]; 
  float val_pfiso_photon03 = ll->pho_pfiso_myphoton03[photon_index];                     
  float val_pfiso_neutral03 = ll->pho_pfiso_myneutral03[photon_index];                   
  float val_pfiso_charged04 = (*ll->pho_pfiso_mycharged04)[photon_index][vtxind];      
  float val_pfiso_photon04 = ll->pho_pfiso_myphoton04[photon_index];                     
  float val_pfiso_neutral04 = ll->pho_pfiso_myneutral04[photon_index];                   
  float val_drtotk_25_99 = ll->pho_drtotk_25_99[photon_index];

  if(MPDEBUGCOPY) cout << "PhotonCiCpfSelectionLevel  02"<<endl;
  float isomax=-99;   int badind=0;
  for(int iv=0; iv<ll->vtx_std_n; iv++) if((*ll->pho_pfiso_mycharged04)[photon_index][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged04)[photon_index][iv]; }
  float val_pfiso_charged_badvtx_04 = (*ll->pho_pfiso_mycharged04)[photon_index][badind];                               
  isomax=-99;  badind=0;
  for(int iv=0; iv<ll->vtx_std_n; iv++) if((*ll->pho_pfiso_mycharged03)[photon_index][iv]>isomax) {badind=iv; isomax=(*ll->pho_pfiso_mycharged03)[photon_index][iv]; }
  float val_pfiso_charged_badvtx_03 = (*ll->pho_pfiso_mycharged03)[photon_index][badind];                               

  if(MPDEBUG) cout<<"DRTOTK: "<<photon_index<<" "<<diphoind<<" "<<ll->pho_drtotk_25_99[photon_index]<<" "
    //MARCO FIX	NOT IMPORTANT	  << ll->dipho_lead_pho_drtotk_25_99[diphoind]<<" "
//MARCO FIX	NOT IMPORTANT		  <<ll->dipho_sublead_pho_drtotk_25_99[diphoind] <<" "
		  <<endl;

  if(MPDEBUGCOPY) cout << "PhotonCiCpfSelectionLevel  03"<<endl;
  float val_sieie = ll->pho_sieie[photon_index];
  float val_hoe = ll->pho_hoe[photon_index];
  float val_r9 = ll->pho_r9[photon_index];
  float val_pixel = (float)ll->pho_haspixseed[photon_index];

  float rhofacpf[6]={0.075, 0.082, 0.143, 0.050, 0.091, 0.106};          //move
  float rhofacbadpf[6]={0.141, 0.149, 0.208, 0.135, 0.162, 0.165};
  float rhofac=rhofacpf[photon_category];
  float val_isosum=val_pfiso_charged03+val_pfiso_photon03-ll->rho*rhofac;
  float rhofacbad=rhofacbadpf[photon_category];
  float val_isosumbad=val_pfiso_charged_badvtx_04+val_pfiso_photon04-ll->rho*rhofacbad;

  if(MPDEBUGCOPY) cout << "PhotonCiCpfSelectionLevel  04"<<endl;
  //cout<<"ind="<<photon_index<<"  vtxind="<<vtxind<<"  badind="<<badind<<"  vtx_std_n="<<ll->vtx_std_n<<"  rho="<<ll->rho<<"  photon_category="<<photon_category<<" rhofac="<<rhofac<<"rhofacbad="<<rhofacbad<<endl;
  if(print) {
    float val_isosumoet=(val_isosum+2.8)*50/val_et;
    float val_isosumoetbad=(val_isosumbad+4.8)*50/val_et;
    cout<<photon_index<<" "<<endl;
    cout<<"val_isosumoet "<<val_isosumoet<<endl;
    cout<<"val_isosumoetbad "<<val_isosumoetbad<<endl;
    cout<<"val_pfiso_charged03 "<<val_pfiso_charged03<<endl;
    cout<<"val_sieie "<<ll->pho_sieie[photon_index]<<endl;
    cout<<"val_hoe "<<ll->pho_hoe[photon_index]<<endl;
    cout<<"val_r9 "<<ll->pho_r9[photon_index]<<endl;
    cout<<"val_drtotk_25_99 "<<ll->pho_drtotk_25_99[photon_index]<<endl;
    
    cout<<"              val_pfiso_charged03 "<<val_pfiso_charged03<<endl;
    cout<<"              val_pfiso_photon03 "<<val_pfiso_photon03<<endl;
    cout<<"              val_pfiso_charged04 "<<val_pfiso_charged04<<endl;
    cout<<"              val_pfiso_charged03 "<<val_pfiso_charged03<<endl;
    
  }
  if(MPDEBUGCOPY) cout << "PhotonCiCpfSelectionLevel  05"<<endl;

  //  Here we use the CiC routines to calculate the photon cut level for either 4 or 6 categories
  //  The parameters are the cut variables so that the choice of the photons, vertex etc. can be separated from the basic Id

  //cout<<"   PhotonCiCSelectionLevel  "<<photon_index<<"  "<<ncat<<"  "<<vtxind<<"  "<<doSublead<<"  "<<diphoind<<"  "<<val_r9<<"  "<<val_eta<<"  "<<val_isosumoet<<"  "
  //    <<val_isosumoetbad<<"  "<<val_trkisooet<<"  "<<val_sieie<<"  "<<val_hoe<<"  "<<val_drtotk_25_99<<endl;

  cutlevelpassed = photonCutLevel6pf(val_r9, val_eta, val_et, val_isosum, val_isosumbad, val_pfiso_charged03, val_pfiso_neutral04, val_sieie, val_hoe, val_drtotk_25_99);
  if(MPDEBUGCOPY) cout << "PhotonCiCpfSelectionLevel  BEGIN"<<endl;


  return cutlevelpassed;

}

