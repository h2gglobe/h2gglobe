#include "../interface/ZMuMuGammaAnalysis.h"

#include "PhotonReducedInfo.h"
#include "Sorters.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "JetAnalysis/interface/JetHandler.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

#define PADEBUG 0

using namespace std;

ZMuMuGammaAnalysis::ZMuMuGammaAnalysis(){}
ZMuMuGammaAnalysis::~ZMuMuGammaAnalysis(){}

void readEnergyScaleOffsetsZMuMuGamma(const std::string &fname, EnergySmearer::energySmearingParameters::eScaleVector &escaleOffsets,
                            EnergySmearer::energySmearingParameters::phoCatVector &photonCategories, bool data=true
                            )
{
    // read in energy scale corrections to be applied in run ranges
    std::fstream in(fname.c_str());
    assert( in );
    char line[200];
    float EBHighR9, EBLowR9, EBm4HighR9, EBm4LowR9, EEHighR9, EELowR9;
    char catname[200];
    float mineta, maxeta, minr9, maxr9, offset, err;
    int type;
    int  first, last;
    do {
        in.getline( line, 200, '\n' );

        if( sscanf(line,"%d %d %f %f %f %f %f %f",&first, &last, &EBHighR9, &EBLowR9, &EBm4HighR9, &EBm4LowR9, &EEHighR9, &EELowR9) == 8 ) {
            std::cerr << "Energy scale by run " <<  first<< " " <<  last<< " " <<  EBHighR9<< " " <<  EBLowR9 << " " <<  EBm4HighR9<< " " <<  EBm4LowR9<< " " <<  EEHighR9<< " " <<  EELowR9 << std::endl;

            assert( ! data );
            escaleOffsets.push_back(EnergyScaleOffset(first,last));
            escaleOffsets.back().scale_offset["EBHighR9"] = -1.*EBHighR9;
            escaleOffsets.back().scale_offset["EBLowR9"]  = -1.*EBLowR9;
            escaleOffsets.back().scale_offset["EBm4HighR9"] = -1.*EBm4HighR9;
            escaleOffsets.back().scale_offset["EBm4LowR9"]  = -1.*EBm4LowR9;
            escaleOffsets.back().scale_offset["EEHighR9"] = -1.*EEHighR9;
            escaleOffsets.back().scale_offset["EELowR9"]  = -1.*EELowR9;
            escaleOffsets.back().scale_offset_error["EBHighR9"] = 0.;
            escaleOffsets.back().scale_offset_error["EBLowR9"]  = 0.;
            escaleOffsets.back().scale_offset_error["EBm4HighR9"] = 0.;
            escaleOffsets.back().scale_offset_error["EBm4LowR9"]  = 0.;
            escaleOffsets.back().scale_offset_error["EEHighR9"] = 0.;
            escaleOffsets.back().scale_offset_error["EELowR9"]  = 0.;
        } else if( sscanf(line,"%s %d %f %f %f %f %d %d %f %f", &catname, &type, &mineta, &maxeta, &minr9, &maxr9, &first, &last, &offset, &err  ) == 10 ) {
        std::cerr << "Energy scale (or smering) by run " <<  catname << " " << type << " " << mineta << " " << maxeta << " " << minr9 << " " << maxr9 << " " << first << " " << last << " " << offset << " " << err << std::endl;

        assert( type>=0 && type<=2 );

            EnergySmearer::energySmearingParameters::eScaleVector::reverse_iterator escaleOffset =
                find(escaleOffsets.rbegin(),escaleOffsets.rend(),std::make_pair(first,last));
            if( escaleOffset == escaleOffsets.rend() ) {
                std::cerr << "  adding new range range " << first << " " << last << std::endl;
                escaleOffsets.push_back(EnergyScaleOffset(first,last));
                escaleOffset = escaleOffsets.rbegin();
            }
            // chck if the category is already defined
            if( find(photonCategories.begin(), photonCategories.end(), std::string(catname) ) == photonCategories.end() ) {
                std::cerr << "  defining new category" << std::endl;
                photonCategories.push_back(PhotonCategory(mineta,maxeta,minr9,maxr9,(PhotonCategory::photon_type_t)type,catname));
            }
            // assign the scale offset and error for this category and this run range
            escaleOffset->scale_offset[catname] = data ? -offset : offset;
            escaleOffset->scale_offset_error[catname] = err;
        }

    } while( in );

    in.close();
}

// ----------------------------------------------------------------------------------------------------
void ZMuMuGammaAnalysis::Init(LoopAll& l)
{
  if(PADEBUG)
    cout << "InitRealZMuMuAnalysis START"<<endl;

  if(energyCorrectionMethod=="DaunceyAndKenzie"){
    energyCorrected     = (l.pho_residCorrEnergy);
    energyCorrectedError= (l.pho_residCorrResn);
  }else if(energyCorrectionMethod=="Bendavid"){
    energyCorrected     = (l.pho_regr_energy);
    energyCorrectedError= (l.pho_regr_energyerr);
  }else if(energyCorrectionMethod=="BendavidOTF"){
    energyCorrected     = (l.pho_regr_energy_otf);
    energyCorrectedError= (l.pho_regr_energyerr_otf);

    //  }else if(energyCorrectionMethod=="PFRegression"){
  }else{
    assert(doEcorrectionSmear==false);
  }
  if (doEcorrectionSmear) std::cout << "using energy correction type: " << energyCorrectionMethod << std::endl;
  else                    std::cout << "NOT using energy correction (sbattogiu)"<< std::endl;

  if( vtxVarNames.empty() ) {
    vtxVarNames.push_back("ptbal"), vtxVarNames.push_back("ptasym"), vtxVarNames.push_back("logsumpt2");
  }

  /// // trigger

  // /cdaq/physics/Run2011/5e32/v4.2/HLT/V2
  triggerSelections.push_back(TriggerSelection(160404,161176));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");

  // /cdaq/physics/Run2011/5e32/v6.1/HLT/V1
  triggerSelections.push_back(TriggerSelection(161216,165633));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon20_R9Id_Photon18_R9Id_v");

  // /cdaq/physics/Run2011/1e33/v2.3/HLT/V1
  triggerSelections.push_back(TriggerSelection(165970,166967));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

  // /cdaq/physics/Run2011/1e33/v2.3/HLT/V1
  triggerSelections.push_back(TriggerSelection(167039,173198));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

  // /cdaq/physics/Run2011/1e33/v2.3/HLT/V1
  triggerSelections.push_back(TriggerSelection(165970,166967));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

  // /cdaq/physics/Run2011/1e33/v2.3/HLT/V3
  triggerSelections.push_back(TriggerSelection(167039,173198));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

  // /cdaq/physics/Run2011/3e33/v1.1/HLT/V1
  triggerSelections.push_back(TriggerSelection(173236,178380));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

  // /cdaq/physics/Run2011/5e33/v1.4/HLT/V3
  triggerSelections.push_back(TriggerSelection(178420,190455));
  triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

  triggerSelections.push_back(TriggerSelection(190456,194269));
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");

  triggerSelections.back().addpath("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v");

  triggerSelections.back().addpath("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_Photon22_R9Id85_v");

  triggerSelections.push_back(TriggerSelection(194270,-1));
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");

  triggerSelections.back().addpath("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v");
  triggerSelections.back().addpath("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v");

  triggerSelections.back().addpath("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  triggerSelections.back().addpath("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
  triggerSelections.back().addpath("HLT_Photon36_R9Id85_Photon22_R9Id85_v");

  // n-1 plots for VBF tag 2011
  l.SetCutVariables("cut_VBFLeadJPt",         &myVBFLeadJPt);
  l.SetCutVariables("cut_VBFSubJPt",          &myVBFSubJPt);
  l.SetCutVariables("cut_VBF_Mjj",            &myVBF_Mjj);
  l.SetCutVariables("cut_VBF_dEta",           &myVBFdEta);
  l.SetCutVariables("cut_VBF_Zep",            &myVBFZep);
  l.SetCutVariables("cut_VBF_dPhi",           &myVBFdPhi);
  l.SetCutVariables("cut_VBF_Mgg0",           &myVBF_Mgg);
  l.SetCutVariables("cut_VBF_Mgg2",           &myVBF_Mgg);
  l.SetCutVariables("cut_VBF_Mgg4",           &myVBF_Mgg);
  l.SetCutVariables("cut_VBF_Mgg10",          &myVBF_Mgg);
  l.SetCutVariables("cut_VBF_Mgg4_100_180",   &myVBF_Mgg);
  l.SetCutVariables("cut_VBF_Mgg2_100_180",   &myVBF_Mgg);

  if( vbfVsDiphoVbfSelection ) {
	multiclassVbfSelection = true;
	assert(mvaVbfCatBoundaries.empty() );
	mvaVbfCatBoundaries = multiclassVbfCatBoundaries0;
  }
  if( mvaVbfSelection || multiclassVbfSelection || bookDiPhoCutsInVbf  ) {
    l.SetCutVariables("cut_VBF_DiPhoPtOverM",   &myVBFDiPhoPtOverM);
    l.SetCutVariables("cut_VBF_LeadPhoPtOverM", &myVBFLeadPhoPtOverM);
    l.SetCutVariables("cut_VBF_SubPhoPtOverM",  &myVBFSubPhoPtOverM);
  }

  if( mvaVbfSelection || multiclassVbfSelection ) {

	tmvaVbfReader_ = new TMVA::Reader( "!Color:!Silent" );

	tmvaVbfReader_->AddVariable("jet1pt"              , &myVBFLeadJPt);
	tmvaVbfReader_->AddVariable("jet2pt"	          , &myVBFSubJPt);
	tmvaVbfReader_->AddVariable("abs(jet1eta-jet2eta)", &myVBFdEta);
	tmvaVbfReader_->AddVariable("mj1j2"		  , &myVBF_Mjj);
	tmvaVbfReader_->AddVariable("zepp"		  , &myVBFZep);
	tmvaVbfReader_->AddVariable("dphi"		  , &myVBFdPhi);
	if( mvaVbfUseDiPhoPt ) {
      tmvaVbfReader_->AddVariable("diphopt/diphoM"      , &myVBFDiPhoPtOverM);
	}
	if( mvaVbfUsePhoPt   ) {
      tmvaVbfReader_->AddVariable("pho1pt/diphoM"	  , &myVBFLeadPhoPtOverM);
      tmvaVbfReader_->AddVariable("pho2pt/diphoM"       , &myVBFSubPhoPtOverM);
	}

	tmvaVbfReader_->BookMVA( mvaVbfMethod, mvaVbfWeights );

  }




  // n-1 plots for VH hadronic tag 2011
  l.SetCutVariables("cut_VHhadLeadJPt",      &myVHhadLeadJPt);
  l.SetCutVariables("cut_VHhadSubJPt",       &myVHhadSubJPt);
  l.SetCutVariables("cut_VHhad_Mjj",         &myVHhad_Mjj);
  l.SetCutVariables("cut_VHhad_dEta",        &myVHhaddEta);
  l.SetCutVariables("cut_VHhad_Zep",         &myVHhadZep);
  l.SetCutVariables("cut_VHhad_dPhi",        &myVHhaddPhi);
  l.SetCutVariables("cut_VHhad_Mgg0",        &myVHhad_Mgg);
  l.SetCutVariables("cut_VHhad_Mgg2",        &myVHhad_Mgg);
  l.SetCutVariables("cut_VHhad_Mgg4",        &myVHhad_Mgg);
  l.SetCutVariables("cut_VHhad_Mgg10",        &myVHhad_Mgg);
  l.SetCutVariables("cut_VHhad_Mgg2_100_160",        &myVHhad_Mgg);
  l.SetCutVariables("cut_VHhad_Mgg4_100_160",        &myVHhad_Mgg);

  // n-1 plot for ClassicCats
  l.SetCutVariables("cutnm1hir9EB_r9",             &sublead_r9);
  l.SetCutVariables("cutnm1hir9EB_isoOverEt",      &sublead_isoOverEt);
  l.SetCutVariables("cutnm1hir9EB_badisoOverEt",   &sublead_badisoOverEt);
  l.SetCutVariables("cutnm1hir9EB_trkisooet",      &sublead_trkisooet);
  l.SetCutVariables("cutnm1hir9EB_sieie",          &sublead_sieie);
  l.SetCutVariables("cutnm1hir9EB_drtotk",         &sublead_drtotk);
  l.SetCutVariables("cutnm1hir9EB_hovere",         &sublead_hovere);
  l.SetCutVariables("cutnm1hir9EB_Mgg",            &sublead_mgg);

  l.SetCutVariables("cutnm1lor9EB_r9",             &sublead_r9);
  l.SetCutVariables("cutnm1lor9EB_isoOverEt",      &sublead_isoOverEt);
  l.SetCutVariables("cutnm1lor9EB_badisoOverEt",   &sublead_badisoOverEt);
  l.SetCutVariables("cutnm1lor9EB_trkisooet",      &sublead_trkisooet);
  l.SetCutVariables("cutnm1lor9EB_sieie",          &sublead_sieie);
  l.SetCutVariables("cutnm1lor9EB_drtotk",         &sublead_drtotk);
  l.SetCutVariables("cutnm1lor9EB_hovere",         &sublead_hovere);
  l.SetCutVariables("cutnm1lor9EB_Mgg",            &sublead_mgg);

  l.SetCutVariables("cutnm1hir9EE_r9",             &sublead_r9);
  l.SetCutVariables("cutnm1hir9EE_isoOverEt",      &sublead_isoOverEt);
  l.SetCutVariables("cutnm1hir9EE_badisoOverEt",   &sublead_badisoOverEt);
  l.SetCutVariables("cutnm1hir9EE_trkisooet",      &sublead_trkisooet);
  l.SetCutVariables("cutnm1hir9EE_sieie",          &sublead_sieie);
  l.SetCutVariables("cutnm1hir9EE_drtotk",         &sublead_drtotk);
  l.SetCutVariables("cutnm1hir9EE_hovere",         &sublead_hovere);
  l.SetCutVariables("cutnm1hir9EE_Mgg",            &sublead_mgg);

  l.SetCutVariables("cutnm1lor9EE_r9",             &sublead_r9);
  l.SetCutVariables("cutnm1lor9EE_isoOverEt",      &sublead_isoOverEt);
  l.SetCutVariables("cutnm1lor9EE_badisoOverEt",   &sublead_badisoOverEt);
  l.SetCutVariables("cutnm1lor9EE_trkisooet",      &sublead_trkisooet);
  l.SetCutVariables("cutnm1lor9EE_sieie",          &sublead_sieie);
  l.SetCutVariables("cutnm1lor9EE_drtotk",         &sublead_drtotk);
  l.SetCutVariables("cutnm1lor9EE_hovere",         &sublead_hovere);
  l.SetCutVariables("cutnm1lor9EE_Mgg",            &sublead_mgg);

  if(includeVHlep) {
    l.SetCutVariables("cutEl_leptonSig",    &myEl_leptonSig);
    l.SetCutVariables("cutEl_elpt",         &myEl_elpt);
    l.SetCutVariables("cutEl_oEsuboP",      &myEl_oEsuboP);
    l.SetCutVariables("cutEl_D0",           &myEl_D0     );
    l.SetCutVariables("cutEl_DZ",           &myEl_DZ     );
    l.SetCutVariables("cutEl_mishit",       &myEl_mishit );
    l.SetCutVariables("cutEl_conv",         &myEl_conv   );
    l.SetCutVariables("cutEl_detain",       &myEl_detain );
    l.SetCutVariables("cutEl_dphiin",       &myEl_dphiin );
    l.SetCutVariables("cutEl_sieie",        &myEl_sieie  );
    l.SetCutVariables("cutEl_sieie2",       &myEl_sieie  );
    l.SetCutVariables("cutEl_hoe",          &myEl_hoe    );
    l.SetCutVariables("cutEl_drlead",       &myEl_drlead );
    l.SetCutVariables("cutEl_drsub",        &myEl_drsub  );
    l.SetCutVariables("cutEl_melead",       &myEl_melead );
    l.SetCutVariables("cutEl_meleadveto10", &myEl_meleadveto10 );
    l.SetCutVariables("cutEl_meleadveto15", &myEl_meleadveto15 );
    l.SetCutVariables("cutEl_mesub",        &myEl_mesub  );
    l.SetCutVariables("cutEl_mesubveto5",   &myEl_mesubveto5  );
    l.SetCutVariables("cutEl_mesubveto10",  &myEl_mesubveto10  );
    l.SetCutVariables("cutEl_reliso",       &myEl_reliso );
    l.SetCutVariables("cutEl_iso",          &myEl_iso    );
    l.SetCutVariables("cutEl_mvaNonTrig",   &myEl_mvaNonTrig);
    l.SetCutVariables("cutEl_dZ_ee",        &myEl_dZ_ee);
    l.SetCutVariables("cutEl_mass_ee",      &myEl_mass_ee);
    l.SetCutVariables("cutEl_inwindow_ee",  &myEl_inwindow_ee);
    l.SetCutVariables("cutEl_ptlead",       &myEl_ptlead    );
    l.SetCutVariables("cutEl_ptsub",        &myEl_ptsub     );
    l.SetCutVariables("cutEl_ptleadom",       &myEl_ptleadom    );
    l.SetCutVariables("cutEl_ptsubom",        &myEl_ptsubom     );
    l.SetCutVariables("cutEl_elvetolead",   &myEl_elvetolead);
    l.SetCutVariables("cutEl_elvetosub",    &myEl_elvetosub );
    l.SetCutVariables("cutEl_ptgg",         &myEl_ptgg      );
    l.SetCutVariables("cutEl_phomaxeta",    &myEl_phomaxeta );
    l.SetCutVariables("cutEl_sumpt3",       &myEl_sumpt3    );
    l.SetCutVariables("cutEl_sumpt4",       &myEl_sumpt4    );
    l.SetCutVariables("cutEl_dRtklead",     &myEl_dRtklead  );
    l.SetCutVariables("cutEl_dRtksub",      &myEl_dRtksub   );
    l.SetCutVariables("cutEl_MVAlead",      &myEl_MVAlead   );
    l.SetCutVariables("cutEl_MVAsub",       &myEl_MVAsub    );
    l.SetCutVariables("cutEl_diphomva",     &myEl_diphomva  );
    l.SetCutVariables("cutEl_CiClead",      &myEl_CiClead   );
    l.SetCutVariables("cutEl_CiCsub",       &myEl_CiCsub    );
    l.SetCutVariables("cutEl_mgg",          &myEl_mgg       );
    l.SetCutVariables("cutEl_MET",          &myEl_MET       );
    l.SetCutVariables("cutEl_METphi",       &myEl_METphi    );
    l.SetCutVariables("cutEl_presellead",   &myEl_presellead );
    l.SetCutVariables("cutEl_matchellead",  &myEl_matchellead);
    l.SetCutVariables("cutEl_preselsub",    &myEl_preselsub  );
    l.SetCutVariables("cutEl_matchelsub",   &myEl_matchelsub );
    l.SetCutVariables("cutEl_category",     &myEl_category );
    l.SetCutVariables("cutEl_ElePho",       &myEl_ElePho );
    l.SetCutVariables("cutEl_passelcuts",   &myEl_passelcuts );
  }


  // CiC initialization
  // FIXME should move this to GeneralFunctions
  l.runCiC = true;
  const int phoNCUTS = LoopAll::phoNCUTS;
  const int phoCiC6NCATEGORIES = LoopAll::phoCiC6NCATEGORIES;
  const int phoCiC4NCATEGORIES = LoopAll::phoCiC4NCATEGORIES;
  const int phoNCUTLEVELS = LoopAll::phoNCUTLEVELS;

  for(int iLevel=0; iLevel<phoNCUTLEVELS; ++iLevel) {
    float cic6_cuts_lead[phoNCUTS][phoCiC6NCATEGORIES];
    float cic6_cuts_sublead[phoNCUTS][phoCiC6NCATEGORIES];
    float cic4_cuts_lead[phoNCUTS][phoCiC4NCATEGORIES];
    float cic4_cuts_sublead[phoNCUTS][phoCiC4NCATEGORIES];
    float cic4pf_cuts_lead[phoNCUTS][phoCiC4NCATEGORIES];
    float cic4pf_cuts_sublead[phoNCUTS][phoCiC4NCATEGORIES];

    // get the cut values for the current photon CIC level
    l.SetPhotonCutsInCategories((LoopAll::phoCiCIDLevel)iLevel,
                                &cic6_cuts_lead[0][0], &cic6_cuts_sublead[0][0],
                                &cic4_cuts_lead[0][0], &cic4_cuts_sublead[0][0],
                                &cic4pf_cuts_lead[0][0], &cic4pf_cuts_sublead[0][0]);

    // rearrange the returned values to arrays with more meaningful names
    float * cic6_cuts_arrays_lead[phoNCUTS] = {
      &l.cic6_cut_lead_isosumoet[0][0],
      &l.cic6_cut_lead_isosumoetbad[0][0],
      &l.cic6_cut_lead_trkisooet[0][0],
      &l.cic6_cut_lead_sieie[0][0],
      &l.cic6_cut_lead_hovere[0][0],
      &l.cic6_cut_lead_r9[0][0],
      &l.cic6_cut_lead_drtotk_25_99[0][0],
      &l.cic6_cut_lead_pixel[0][0]
    };

    float * cic6_cuts_arrays_sublead[phoNCUTS] = {
      &l.cic6_cut_sublead_isosumoet[0][0],
      &l.cic6_cut_sublead_isosumoetbad[0][0],
      &l.cic6_cut_sublead_trkisooet[0][0],
      &l.cic6_cut_sublead_sieie[0][0],
      &l.cic6_cut_sublead_hovere[0][0],
      &l.cic6_cut_sublead_r9[0][0],
      &l.cic6_cut_sublead_drtotk_25_99[0][0],
      &l.cic6_cut_sublead_pixel[0][0]
    };

    float * cic4_cuts_arrays_lead[phoNCUTS] = {
      &l.cic4_cut_lead_isosumoet[0][0],
      &l.cic4_cut_lead_isosumoetbad[0][0],
      &l.cic4_cut_lead_trkisooet[0][0],
      &l.cic4_cut_lead_sieie[0][0],
      &l.cic4_cut_lead_hovere[0][0],
      &l.cic4_cut_lead_r9[0][0],
      &l.cic4_cut_lead_drtotk_25_99[0][0],
      &l.cic4_cut_lead_pixel[0][0]
    };

    float * cic4_cuts_arrays_sublead[phoNCUTS] = {
      &l.cic4_cut_sublead_isosumoet[0][0],
      &l.cic4_cut_sublead_isosumoetbad[0][0],
      &l.cic4_cut_sublead_trkisooet[0][0],
      &l.cic4_cut_sublead_sieie[0][0],
      &l.cic4_cut_sublead_hovere[0][0],
      &l.cic4_cut_sublead_r9[0][0],
      &l.cic4_cut_sublead_drtotk_25_99[0][0],
      &l.cic4_cut_sublead_pixel[0][0]
    };

    float * cic4pf_cuts_arrays_lead[phoNCUTS] = {
      &l.cic4pf_cut_lead_isosumoet[0][0],
      &l.cic4pf_cut_lead_isosumoetbad[0][0],
      &l.cic4pf_cut_lead_trkisooet[0][0],
      &l.cic4pf_cut_lead_sieie[0][0],
      &l.cic4pf_cut_lead_hovere[0][0],
      &l.cic4pf_cut_lead_r9[0][0],
      &l.cic4pf_cut_lead_drtotk_25_99[0][0],
      &l.cic4pf_cut_lead_pixel[0][0]
    };

    float * cic4pf_cuts_arrays_sublead[phoNCUTS] = {
      &l.cic4pf_cut_sublead_isosumoet[0][0],
      &l.cic4pf_cut_sublead_isosumoetbad[0][0],
      &l.cic4pf_cut_sublead_trkisooet[0][0],
      &l.cic4pf_cut_sublead_sieie[0][0],
      &l.cic4pf_cut_sublead_hovere[0][0],
      &l.cic4pf_cut_sublead_r9[0][0],
      &l.cic4pf_cut_sublead_drtotk_25_99[0][0],
      &l.cic4pf_cut_sublead_pixel[0][0]
    };
    for(int iCut=0; iCut<phoNCUTS; ++iCut) {
      for(int iCat=0; iCat<phoCiC6NCATEGORIES; ++iCat) {
        cic6_cuts_arrays_lead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_lead[iCut][iCat];
        cic6_cuts_arrays_sublead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_sublead[iCut][iCat];
      }
      for(int iCat=0; iCat<phoCiC4NCATEGORIES; ++iCat) {
        cic4_cuts_arrays_lead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_lead[iCut][iCat];
        cic4_cuts_arrays_sublead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_sublead[iCut][iCat];
      }
      for(int iCat=0; iCat<phoCiC4NCATEGORIES; ++iCat) {
        cic4pf_cuts_arrays_lead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4pf_cuts_lead[iCut][iCat];
        cic4pf_cuts_arrays_sublead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4pf_cuts_sublead[iCut][iCat];
      }
    }
  } // end of loop over all photon cut levels

    //--------------------

  if( tmvaPerVtxWeights != ""  ) {
	if( tmvaPerVtxVariables.empty() ) {
      tmvaPerVtxVariables.push_back("ptbal"), tmvaPerVtxVariables.push_back("ptasym"), tmvaPerVtxVariables.push_back("logsumpt2");
      if( addConversionToMva ) {
		tmvaPerVtxVariables.push_back("limPullToConv");
		tmvaPerVtxVariables.push_back("nConv");
      }
	}
    tmvaPerVtxReader_ = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookVariables( *tmvaPerVtxReader_, tmvaPerVtxVariables );
    tmvaPerVtxReader_->BookMVA( tmvaPerVtxMethod, tmvaPerVtxWeights );
  } else {
    tmvaPerVtxReader_ = 0;
  }
  if( tmvaPerEvtWeights != "" ) {
    tmvaPerEvtReader_ = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookPerEventVariables( *tmvaPerEvtReader_ );
    tmvaPerEvtReader_->BookMVA( tmvaPerEvtMethod, tmvaPerEvtWeights );
  } else {
    tmvaPerEvtReader_ = 0;
  }
  assert( !mvaVertexSelection || tmvaPerVtxReader_ != 0 );

  eSmearDataPars.categoryType = "2CatR9_EBEBm4EE";
  eSmearDataPars.byRun = true;
  //eSmearDataPars.n_categories = 4; //GF
  eSmearDataPars.n_categories = 6; //GF
  std::cerr << "Reading energy scale offsets " << scale_offset_file << std::endl;
  readEnergyScaleOffsetsZMuMuGamma(scale_offset_file, eSmearDataPars.scale_offset_byrun, eSmearDataPars.photon_categories);
  // if the scale offset file defines the categories set the category type to automatic
  if( ! eSmearDataPars.photon_categories.empty() ) {
    eSmearDataPars.categoryType = "Automagic";
    eSmearDataPars.n_categories = -1;
  }
  // E resolution smearing NOT applied to data
  eSmearDataPars.smearing_sigma["EBHighR9"] = 0.;
  eSmearDataPars.smearing_sigma["EBLowR9"]  = 0.;
  eSmearDataPars.smearing_sigma["EBm4HighR9"] = 0.;
  eSmearDataPars.smearing_sigma["EBm4LowR9"]  = 0.;
  eSmearDataPars.smearing_sigma["EEHighR9"] = 0.;
  eSmearDataPars.smearing_sigma["EELowR9"]  = 0.;
  // E resolution systematics NOT applied to data
  eSmearDataPars.smearing_sigma_error["EBHighR9"] = 0.;
  eSmearDataPars.smearing_sigma_error["EBLowR9"]  = 0.;
  eSmearDataPars.smearing_sigma_error["EBm4HighR9"] = 0.;
  eSmearDataPars.smearing_sigma_error["EBm4LowR9"]  = 0.;
  eSmearDataPars.smearing_sigma_error["EEHighR9"] = 0.;
  eSmearDataPars.smearing_sigma_error["EELowR9"]  = 0.;

  // energy scale corrections to Data
  eScaleDataSmearer = new EnergySmearer( eSmearDataPars );
  eScaleDataSmearer->name("E_scale_data");
  eScaleDataSmearer->doEnergy(true);
  eScaleDataSmearer->scaleOrSmear(true);

  if( scale_offset_error_file.empty() ) {
    //eSmearPars.categoryType = "2CatR9_EBEE"; //GF
    eSmearPars.categoryType = "2CatR9_EBEBm4EE";
    eSmearPars.byRun = false;
    //eSmearPars.n_categories = 4; //GF
    eSmearPars.n_categories = 6;
    // E scale is shifted for data, NOT for MC
    eSmearPars.scale_offset["EBHighR9"] = 0.;
    eSmearPars.scale_offset["EBLowR9"]  = 0.;
    eSmearPars.scale_offset["EBm4HighR9"] = 0.;
    eSmearPars.scale_offset["EBm4LowR9"]  = 0.;
    eSmearPars.scale_offset["EEHighR9"] = 0.;
    eSmearPars.scale_offset["EELowR9"]  = 0.;
    // E scale systematics are applied to MC, NOT to data
    eSmearPars.scale_offset_error["EBHighR9"] = scale_offset_error_EBHighR9;
    eSmearPars.scale_offset_error["EBLowR9"]  = scale_offset_error_EBLowR9;
    eSmearPars.scale_offset_error["EBm4HighR9"] = scale_offset_error_EBHighR9;
    eSmearPars.scale_offset_error["EBm4LowR9"]  = scale_offset_error_EBLowR9;
    eSmearPars.scale_offset_error["EEHighR9"] = scale_offset_error_EEHighR9;
    eSmearPars.scale_offset_error["EELowR9"]  = scale_offset_error_EELowR9;
    // E resolution smearing applied to MC
    eSmearPars.smearing_sigma["EBHighR9"] = smearing_sigma_EBHighR9;
    eSmearPars.smearing_sigma["EBLowR9"]  = smearing_sigma_EBLowR9;
    eSmearPars.smearing_sigma["EBm4HighR9"] = smearing_sigma_EBm4HighR9;
    eSmearPars.smearing_sigma["EBm4LowR9"]  = smearing_sigma_EBm4LowR9;
    eSmearPars.smearing_sigma["EEHighR9"] = smearing_sigma_EEHighR9;
    eSmearPars.smearing_sigma["EELowR9"]  = smearing_sigma_EELowR9;
    // E resolution systematics applied to MC
    eSmearPars.smearing_sigma_error["EBHighR9"] = smearing_sigma_error_EBHighR9;
    eSmearPars.smearing_sigma_error["EBLowR9"]  = smearing_sigma_error_EBLowR9;
    eSmearPars.smearing_sigma_error["EBm4HighR9"] = smearing_sigma_error_EBm4HighR9;
    eSmearPars.smearing_sigma_error["EBm4LowR9"]  = smearing_sigma_error_EBm4LowR9;
    eSmearPars.smearing_sigma_error["EEHighR9"] = smearing_sigma_error_EEHighR9;
    eSmearPars.smearing_sigma_error["EELowR9"]  = smearing_sigma_error_EELowR9;
    // error on photon corrections set to a fraction of the correction itself; number below is tentative (GF: push it to .dat)
    eSmearPars.corrRelErr  = 0.5;
  } else {
    // Read energy scale errors and energy smaerings from dat files
    assert( ! scale_offset_error_file.empty() && ! smearing_file.empty() );

    // Use the same format used for the run-dependent energy corrections
    EnergySmearer::energySmearingParameters::eScaleVector tmp_scale_offset, tmp_smearing;
    EnergySmearer::energySmearingParameters::phoCatVector tmp_scale_cat, tmp_smearing_cat;
    readEnergyScaleOffsetsZMuMuGamma(scale_offset_error_file, tmp_scale_offset, tmp_scale_cat,false);
    readEnergyScaleOffsetsZMuMuGamma(smearing_file, tmp_smearing, tmp_smearing_cat,false);

    // make sure that the scale correction and smearing info is as expected
    assert( tmp_scale_offset.size() == 1); assert( tmp_smearing.size() == 1 );
    assert( ! tmp_smearing_cat.empty() );
    /// assert( tmp_smearing_cat == tmp_scale_cat );

    // copy the read info to the smarer parameters
    eSmearPars.categoryType = "Automagic";
    eSmearPars.byRun = false;
    eSmearPars.n_categories = tmp_smearing_cat.size();
    eSmearPars.photon_categories = tmp_smearing_cat;

    eSmearPars.scale_offset = tmp_scale_offset[0].scale_offset;
    eSmearPars.scale_offset_error = tmp_scale_offset[0].scale_offset_error;

    eSmearPars.smearing_sigma = tmp_smearing[0].scale_offset;
    eSmearPars.smearing_sigma_error = tmp_smearing[0].scale_offset_error;
  }

  // energy scale systematics to MC
  eScaleSmearer = new EnergySmearer( eSmearPars );
  eScaleSmearer->name("E_scale");
  eScaleSmearer->doEnergy(true);
  eScaleSmearer->scaleOrSmear(true);
  eScaleSmearer->syst_only(true);

  eResolSmearer = new EnergySmearer( eSmearPars );
  eResolSmearer->name("E_res");
  eResolSmearer->doEnergy(false);
  eResolSmearer->scaleOrSmear(false);

  if( doEcorrectionSmear ) {
    eCorrSmearer = new EnergySmearer( eSmearPars );
    eCorrSmearer->name("E_corr");
    // activating pho corrections to this instance of EnergySmearer, implies that it won't touch Escale and Eresolution
    eCorrSmearer->doCorrections(true);
  }

  if (l.typerun == 2 || l.typerun == 1) {
  }
  // MassResolution
  massResolutionCalculator = new MassResolution();
  /* -------------------------------------------------------------------------------------------
     Pileup Reweighting
     https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
     ----------------------------------------------------------------------------------------------  */
  if (puHist != "" && puHist != "auto" ) {
    if(PADEBUG)
      cout << "Opening PU file"<<endl;
    TFile* puFile = TFile::Open( puHist );
    if (puFile) {
      TH1 * target = 0;

      if( puTarget != "" ) {
        TFile * puTargetFile = TFile::Open( puTarget );
        assert( puTargetFile != 0 );
        target = (TH1*)puTargetFile->Get("pileup");
        if( target == 0 ) { target = (TH1*)puTargetFile->Get("target_pileup"); }
        target->Scale( 1. / target->Integral() );
      }

      if( puMap != "" ) {
        loadPuMap(puMap, puFile, target);
      } else {
        loadPuWeights(0, puFile, target);
      }
      puFile->Close();
    }
    else {
      cout<<"Error opening " <<puHist<<" pileup reweighting histogram, using 1.0"<<endl;
      weights[0].resize(50);
      for (unsigned int i=0; i<weights[0].size(); i++) weights[0][i] = 1.0;
    }
    if(PADEBUG)
      cout << "Opening PU file END"<<endl;
  } else if ( puHist == "auto" ) {
	TFile * puTargetFile = TFile::Open( puTarget );
	assert( puTargetFile != 0 );
	puTargetHist = (TH1*)puTargetFile->Get("pileup");
	if( puTargetHist == 0 ) { puTargetHist = (TH1*)puTargetFile->Get("target_pileup"); }
	puTargetHist = (TH1*)puTargetHist->Clone();
	puTargetHist->SetDirectory(0);
	puTargetHist->Scale( 1. / puTargetHist->Integral() );
	puTargetFile->Close();
  }

  if( recomputeBetas || recorrectJets || rerunJetMva || recomputeJetWp || applyJer || applyJecUnc || l.typerun != l.kFill ) {
	std::cout << "JetHandler: \n"
              << "recomputeBetas " << recomputeBetas << "\n"
              << "recorrectJets " << recorrectJets << "\n"
              << "rerunJetMva " << rerunJetMva << "\n"
              << "recomputeJetWp " << recomputeJetWp
              << std::endl;
	jetHandler_ = new JetHandler(jetHandlerCfg, l);
  }

  if( emulateBeamspot || reweighBeamspot ) {
	assert( emulatedBeamspotWidth != 0. );
	beamspotWidth = emulatedBeamspotWidth;
  }
  if( beamspotWidth == 0. ) {
	beamspotWidth = (dataIs2011 ? 5.8 : 4.8);
  }

  // Load up instances of PhotonFix for local coordinate calculations
  /*
    PhotonFix::initialise("4_2",photonFixDat);
    std::cout << "Regression corrections from -> " << regressionFile.c_str() << std::endl;
    fgbr = new TFile(regressionFile.c_str(),"READ");
    fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
    fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");
    fReaderee = (GBRForest*)fgbr->Get("EECorrection");
    fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");
    fgbr->Close();
  */

  l.SetAllMVA();
  cout << "Weights file is: " << photonLevelNewIDMVA_EB.c_str() << endl;
  l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevelNewIDMVA_EB.c_str());
  l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevelNewIDMVA_EE.c_str());
  
  // --------------------------------------------------------------------
  if(PADEBUG)
    cout << "InitRealZMuMuGammaAnalysis END"<<endl;

  // FIXME book of additional variables

}

//----------------------------------------------------------------------------------------------------
bool ZMuMuGammaAnalysis::SkimEvents(LoopAll& l, int jentry)
{
    if( dataIs2011 ) { l.version=12; }
    return true;
}

// ----------------------------------------------------------------------------------------------------
bool ZMuMuGammaAnalysis::SelectEventsReduction(LoopAll& l, int jentry)
{

    if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
    // require at least two reconstructed photons to store the event

    //if( pho_acc.size() < 2 ) { return false; }

    vtxAna_.clear();
    l.vtx_std_ranked_list->clear();
    l.dipho_vtx_std_sel->clear();
    l.vtx_std_ranked_list->clear();
    l.vtx_std_evt_mva->clear();
    l.vtx_std_sel=0;
    l.pho_mitmva->clear();
    float maxSumPt = 0.;
    l.dipho_n = 0;
    bool oneKinSelected = true;

    // fill ID variables
    if( forcedRho >= 0. ) {
        l.rho = forcedRho;
    } else if ( l.rho == 0. ) {
        l.rho = l.rho_algo1;
    }
    l.FillCICInputs();
    if(reComputeCiCPF) { l.FillCICPFInputs(); }
    l.FillCIC();
    l.FillMuonGsfTracks();

    //Calculate cluster shape variables prior to shape rescaling
    for (int ipho=0;ipho<l.pho_n;ipho++){
      l.pho_s4ratio[ipho] = l.pho_e2x2[ipho]/l.bc_s25[l.sc_bcseedind[l.pho_scind[ipho]]];
      float rr2=l.pho_eseffsixix[ipho]*l.pho_eseffsixix[ipho]+l.pho_eseffsiyiy[ipho]*l.pho_eseffsiyiy[ipho];
      l.pho_ESEffSigmaRR[ipho] = 0.0;
      if(rr2>0. && rr2<999999.) {
        l.pho_ESEffSigmaRR[ipho] = sqrt(rr2);
      }
    }
    
    for (int ipho=0; ipho<l.pho_n; ipho++) {
      vector<float> MVAValues;
      for(int ivtx=0; ivtx<l.vtx_std_n; ++ivtx) {
        TLorentzVector pho_p4 = l.get_pho_p4(ipho, ivtx);
        MVAValues.push_back(l.photonIDMVANew(ipho, ivtx, pho_p4, "MIT"));
      }
      l.pho_mitmva->push_back(MVAValues);
    }
    
    if(l.itype[l.current]<0) {
        bool foundHiggs=FindHiggsObjects(l);
        if(PADEBUG)  cout << " foundHiggs? "<<foundHiggs<<std::endl;
    } else {
        SetNullHiggs(l);
    }

    if( pho_presel.size() < 2 ) {
        // zero or one photons, can't determine a vertex based on photon pairs
        l.vtx_std_ranked_list->push_back( std::vector<int>() );
        for(int ii=0;ii<l.vtx_std_n; ++ii) { l.vtx_std_ranked_list->back().push_back(ii); }
        l.vtx_std_sel = 0;
    } else {
        // fully combinatorial vertex selection
        std::vector<std::pair<int,int> > diphotons;
        for(size_t ip=0; ip<pho_presel.size(); ++ip) {
            for(size_t jp=ip+1; jp<pho_presel.size(); ++jp) {
                diphotons.push_back( std::make_pair( pho_presel[ip], pho_presel[jp] ) );
            }
        }
        l.dipho_n = 0;
        for(size_t id=0; id<diphotons.size(); ++id ) {

	    if( l.dipho_n >= MAX_DIPHOTONS-1 ) { continue; }
            int ipho1 = diphotons[id].first;
            int ipho2 = diphotons[id].second;

            if(PADEBUG)        cout << " SelectEventsReduction going to fill photon info " << endl;
            PhotonInfo pho1=l.fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions,&corrected_pho_energy[0]);
            PhotonInfo pho2=l.fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions,&corrected_pho_energy[0]);
            if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;

            l.vertexAnalysis(vtxAna_, pho1, pho2 );
            std::vector<int> vtxs = l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames, mvaVertexSelection,
                              tmvaPerVtxReader_, tmvaPerVtxMethod);

            TLorentzVector lead_p4 = l.get_pho_p4( ipho2, vtxs[0], &corrected_pho_energy[0] ).Pt();
            TLorentzVector sublead_p4 = l.get_pho_p4( ipho1, vtxs[0], &corrected_pho_energy[0] ).Pt();

            if(sublead_p4.Pt()  > lead_p4.Pt() ) {
                std::swap( diphotons[id].first,  diphotons[id].second );
                std::swap( lead_p4,  sublead_p4 );
            }

            if( lead_p4.Pt() < presel_scet1 || sublead_p4.Pt() < presel_scet2 ||
                fabs(lead_p4.Eta()) > presel_maxeta || fabs(sublead_p4.Eta()) > presel_maxeta ) {
                vtxAna_.discardLastDipho();
                continue;
            }
	    oneKinSelected = true;

            if( ! l.PhotonMITPreSelection(ipho1, vtxs[0], &corrected_pho_energy[0] )
                || ! l.PhotonMITPreSelection(ipho2, vtxs[0], &corrected_pho_energy[0] ) ) {
                vtxAna_.discardLastDipho();
                continue;
            }

            l.vtx_std_ranked_list->push_back(vtxs);
            if( tmvaPerEvtReader_ ) {
                float vtxEvtMva = vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back() );
                l.vtx_std_evt_mva->push_back(vtxEvtMva);
            }
            if( l.vtx_std_ranked_list->back().size() != 0 && ! useDefaultVertex ) {
                l.dipho_vtx_std_sel->push_back( (l.vtx_std_ranked_list)->back()[0] );
            } else {
                l.dipho_vtx_std_sel->push_back(0);
                std::cerr << "NO VERTEX SELECTED " << l.event << " " << l.run << " " << diphotons[id].first << " " << diphotons[id].second << std::endl;
            }

            l.dipho_leadind[l.dipho_n] = diphotons[id].first;
            l.dipho_subleadind[l.dipho_n] = diphotons[id].second;
            l.dipho_vtxind[l.dipho_n] = l.dipho_vtx_std_sel->back();

            l.dipho_sumpt[l.dipho_n] = lead_p4.Pt() + sublead_p4.Pt();

            if( l.dipho_sumpt[l.dipho_n] > maxSumPt ) {
                l.vtx_std_sel = l.dipho_vtx_std_sel->back();
                maxSumPt = l.dipho_sumpt[l.dipho_n];
            }

            // make sure that vertex analysis indexes are in synch
            assert( l.dipho_n == vtxAna_.pairID(ipho1,ipho2) );

            l.dipho_n++;
        }

       MetCorrections2012( l );
    }

    // Post-process jets and compute beta variables for missing vertexes if needed.
    int highestVtx = ( ! l.dipho_vtx_std_sel->empty() ?
		       *std::max_element(l.dipho_vtx_std_sel->begin(), l.dipho_vtx_std_sel->end()) + 1
		       : 1 );
    for(int ivtx = 0; ivtx<highestVtx; ++ivtx ) {
	postProcessJets(l,ivtx);
    }

    return oneKinSelected;
}

// ----------------------------------------------------------------------------------------------------
void ZMuMuGammaAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree)
{
    if( outputTree ) { vtxAna_.branches(outputTree,"vtx_std_"); }

    l.pho_matchingConv = new  std::vector<int>();
    if( outputTree ) { l.Branch_pho_matchingConv(outputTree); }

    l.vtx_std_evt_mva = new std::vector<float>();
    l.vtx_std_ranked_list = new std::vector<std::vector<int> >();
    l.pho_tkiso_recvtx_030_002_0000_10_01 = new std::vector<std::vector<float> >();
    l.pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01 = new std::vector<std::vector<float> >();
    l.pho_cic6cutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic6passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic6cutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic6passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4cutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4cutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4pfcutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4pfpasscuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4pfcutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4pfpasscuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_mitmva = new std::vector<std::vector<float> >();
    l.dipho_vtx_std_sel =  new std::vector<int>();

    if( outputTree ) {
	//l.Branch_pho_ncrys(outputTree);

	l.Branch_vtx_std_evt_mva(outputTree);
	l.Branch_vtx_std_ranked_list(outputTree);
	l.Branch_vtx_std_sel(outputTree);
	l.Branch_pho_tkiso_recvtx_030_002_0000_10_01(outputTree);
	l.Branch_pho_tkiso_badvtx_040_002_0000_10_01(outputTree);
	l.Branch_pho_tkiso_badvtx_id(outputTree);
	l.Branch_pho_pfiso_charged_badvtx_04(outputTree);
	l.Branch_pho_pfiso_charged_badvtx_id(outputTree);
	l.Branch_pho_ZeeVal_tkiso_recvtx_030_002_0000_10_01(outputTree);
	l.Branch_pho_ZeeVal_tkiso_badvtx_040_002_0000_10_01(outputTree);
	l.Branch_pho_ZeeVal_tkiso_badvtx_id(outputTree);
	l.Branch_pho_mitmva(outputTree);
	l.Branch_pho_drtotk_25_99(outputTree);

	l.Branch_dipho_n(outputTree);
	l.Branch_dipho_leadind(outputTree);
	l.Branch_dipho_subleadind(outputTree);
	l.Branch_dipho_vtxind(outputTree);
	l.Branch_dipho_sumpt(outputTree);

	l.Branch_pho_cic6cutlevel_lead( outputTree );
	l.Branch_pho_cic6passcuts_lead( outputTree );
	l.Branch_pho_cic6cutlevel_sublead( outputTree );
	l.Branch_pho_cic6passcuts_sublead( outputTree );
	l.Branch_pho_cic4cutlevel_lead( outputTree );
	l.Branch_pho_cic4passcuts_lead( outputTree );
	l.Branch_pho_cic4cutlevel_sublead( outputTree );
	l.Branch_pho_cic4passcuts_sublead( outputTree );
	l.Branch_pho_cic4pfcutlevel_lead( outputTree );
	l.Branch_pho_cic4pfpasscuts_lead( outputTree );
	l.Branch_pho_cic4pfcutlevel_sublead( outputTree );
	l.Branch_pho_cic4pfpasscuts_sublead( outputTree );

	l.Branch_pho_genmatched(outputTree);
	l.Branch_pho_regr_energy_otf(outputTree);
	l.Branch_pho_regr_energyerr_otf(outputTree);

	l.Branch_jet_algoPF1_genMatched(outputTree);
	l.Branch_jet_algoPF1_vbfMatched(outputTree);
	l.Branch_jet_algoPF1_genPt(outputTree);
	l.Branch_jet_algoPF1_genDr(outputTree);

	//l.Branch_jet_algoPF2_genMatched(outputTree);
	//l.Branch_jet_algoPF2_vbfMatched(outputTree);
	//l.Branch_jet_algoPF2_genPt(outputTree);
	//l.Branch_jet_algoPF2_genDr(outputTree);

	l.Branch_jet_algoPF3_genMatched(outputTree);
	l.Branch_jet_algoPF3_vbfMatched(outputTree);
	l.Branch_jet_algoPF3_genPt(outputTree);
	l.Branch_jet_algoPF3_genDr(outputTree);

	//correctMETinRED
	l.Branch_shiftMET_pt(outputTree);
	l.Branch_shiftMET_phi(outputTree);
	l.Branch_smearMET_pt(outputTree);
	l.Branch_smearMET_phi(outputTree);
	l.Branch_shiftsmearMET_pt(outputTree);
	l.Branch_shiftsmearMET_phi(outputTree);
	l.Branch_shiftscaleMET_pt(outputTree);
	l.Branch_shiftscaleMET_phi(outputTree);
	l.Branch_shiftMET_eta(outputTree);
	l.Branch_shiftMET_e(outputTree);
	l.Branch_shiftscaleMET_eta(outputTree);
	l.Branch_shiftscaleMET_e(outputTree);
    }

    l.gh_higgs_p4 = new TClonesArray("TLorentzVector", 1);
    l.gh_higgs_p4->Clear();
    ((*l.gh_higgs_p4)[0]) = new TLorentzVector();

    l.gh_pho1_p4 = new TClonesArray("TLorentzVector", 1);
    l.gh_pho1_p4->Clear();
    ((*l.gh_pho1_p4)[0]) = new TLorentzVector();

    l.gh_pho2_p4 = new TClonesArray("TLorentzVector", 1);
    l.gh_pho2_p4->Clear();
    ((*l.gh_pho2_p4)[0]) = new TLorentzVector();

    l.gh_vbfq1_p4 = new TClonesArray("TLorentzVector", 1);
    l.gh_vbfq1_p4->Clear();
    ((*l.gh_vbfq1_p4)[0]) = new TLorentzVector();

    l.gh_vbfq2_p4 = new TClonesArray("TLorentzVector", 1);
    l.gh_vbfq2_p4->Clear();
    ((*l.gh_vbfq2_p4)[0]) = new TLorentzVector();

    l.gh_vh1_p4 = new TClonesArray("TLorentzVector", 1);
    l.gh_vh1_p4->Clear();
    ((*l.gh_vh1_p4)[0]) = new TLorentzVector();

    l.gh_vh2_p4 = new TClonesArray("TLorentzVector", 1);
    l.gh_vh2_p4->Clear();
    ((*l.gh_vh2_p4)[0]) = new TLorentzVector();

//    l.METcorrected = new TClonesArray("TLorentzVector", 1);     //met at analysis step
//    l.METcorrected->Clear();                    //met at analysis step
//    ((*l.METcorrected)[0]) = new TLorentzVector();        //met at analysis step


    if( outputTree ) {
    l.Branch_gh_gen2reco1( outputTree );
    l.Branch_gh_gen2reco2( outputTree );
    l.Branch_gh_vbfq1_pdgid( outputTree );
    l.Branch_gh_vbfq2_pdgid( outputTree );
    l.Branch_gh_vh_pdgid( outputTree );
    l.Branch_gh_vh1_pdgid( outputTree );
    l.Branch_gh_vh2_pdgid( outputTree );
//    l.Branch_METcorrected( outputTree );  //met at analysis step
    l.Branch_gh_higgs_p4( outputTree );
    l.Branch_gh_pho1_p4( outputTree );
    l.Branch_gh_pho2_p4( outputTree );
    l.Branch_gh_vbfq1_p4( outputTree );
    l.Branch_gh_vbfq2_p4( outputTree );
    l.Branch_gh_vh1_p4( outputTree );
    l.Branch_gh_vh2_p4( outputTree );
    }
    l.Branch_mu_glo_hasgsftrack(outputTree);
}
