#include <iostream>
#include <cassert>
#include <string>
#include <fstream>

#include <TFile.h>
#include <TTree.h>

#include <dlfcn.h>

#include "GenericAnalysis.h"
#include "parser.h"
#include "Smearing.h"

#include <vector>

#include <TRint.h>

#include "Plotting.h"

using namespace std;

//----------------------------------------------------------------------
// some constants / parameters
//----------------------------------------------------------------------

/** the name of the tree within the input file */
const string inputTreeName = "opttree";

bool interactive = true;

void readEnergyScaleOffsets(const std::string &fname, EnergySmearer::energySmearingParameters::eScaleVector &escaleOffsets, 
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

//----------------------------------------------------------------------

void bookHistogram(HistoContainer *container, map<string, string> values)
{
  int htyp = getUint(values,"htyp");
  switch (htyp)
    {
    case 0: // 1D histos
      container->Add(getString(values,"name"),
                     getString(values,"xaxis", ""),
                     getString(values,"yaxis", ""),
                     getUint(values,"ncat"),
                     getUint(values,"xbins"),
                     getFloat(values,"xmin"),
                     getFloat(values,"xmax")
                     );
      break;


    case 1: // 2D histos
      container->Add(getString(values,"name"),
                     getString(values,"xaxis", ""),
                     getString(values,"yaxis", ""),
                     getUint(values,"ncat"),

                     getUint(values,"xbins"),
                     getFloat(values,"xmin"),
                     getFloat(values,"xmax"),

                     getUint(values,"ybins"),
                     getFloat(values,"ymin"),
                     getFloat(values,"ymax")
                     );
      break;

      
    default:
      throw std::runtime_error("unknown htyp " + boost::lexical_cast<string>(htyp));
    }
}

void buildBkgModel(RooContainer* rooContainer,  map<string, string> values) {

  std::string postfix = (getUint(values,"dataIs2011")?"":"_8TeV");
  int nCategories_ = getUint(values, "nCategories");
  std::vector<int> bkgPolOrderByCat = getVint(values, "bkgPolOrderByCat");

  // sanity check
  assert(bkgPolOrderByCat.size() == nCategories_);

  rooContainer->AddRealVar("CMS_hgg_pol6_0" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_1" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_2" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_3" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_4" + postfix, -0.01,- 1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol6_5" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_0" + postfix, "@0*@0", "CMS_hgg_pol6_0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_1" + postfix, "@0*@0", "CMS_hgg_pol6_1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_2" + postfix, "@0*@0", "CMS_hgg_pol6_2" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_3" + postfix, "@0*@0", "CMS_hgg_pol6_3" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_4" + postfix, "@0*@0", "CMS_hgg_pol6_4" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol6_5" + postfix, "@0*@0", "CMS_hgg_pol6_4" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_pol5_0" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_1" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_2" + postfix, -0.1, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_3" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddRealVar("CMS_hgg_pol5_4" + postfix, -0.01, -1.0, 1.0);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_0" + postfix, "@0*@0", "CMS_hgg_pol5_0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_1" + postfix, "@0*@0", "CMS_hgg_pol5_1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_2" + postfix, "@0*@0", "CMS_hgg_pol5_2" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_3" + postfix, "@0*@0", "CMS_hgg_pol5_3" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modpol5_4" + postfix, "@0*@0", "CMS_hgg_pol5_4" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_quartic0"+postfix,-0.1,-1.0,1.0);
  rooContainer->AddRealVar("CMS_hgg_quartic1"+postfix,-0.1,-1.0,1.0);
  rooContainer->AddRealVar("CMS_hgg_quartic2"+postfix,-0.1,-1.0,1.0);
  rooContainer->AddRealVar("CMS_hgg_quartic3"+postfix,-0.01,-1.0,1.0);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic0" + postfix, "@0*@0", "CMS_hgg_quartic0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic1" + postfix, "@0*@0", "CMS_hgg_quartic1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic2" + postfix, "@0*@0", "CMS_hgg_quartic2" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquartic3" + postfix, "@0*@0", "CMS_hgg_quartic3" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_quad0" + postfix, -0.1, -1.5, 1.5);
  rooContainer->AddRealVar("CMS_hgg_quad1" + postfix, -0.01, -1.5, 1.5);
  rooContainer->AddFormulaVar("CMS_hgg_modquad0" + postfix, "@0*@0", "CMS_hgg_quad0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modquad1" + postfix, "@0*@0", "CMS_hgg_quad1" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_cubic0" + postfix, -0.1, -1.5, 1.5);
  rooContainer->AddRealVar("CMS_hgg_cubic1" + postfix, -0.1, -1.5, 1.5);
  rooContainer->AddRealVar("CMS_hgg_cubic2" + postfix, -0.01, -1.5, 1.5);
  rooContainer->AddFormulaVar("CMS_hgg_modcubic0" + postfix, "@0*@0", "CMS_hgg_cubic0" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modcubic1" + postfix, "@0*@0", "CMS_hgg_cubic1" + postfix);
  rooContainer->AddFormulaVar("CMS_hgg_modcubic2" + postfix, "@0*@0", "CMS_hgg_cubic2" + postfix);
  
  rooContainer->AddRealVar("CMS_hgg_lin0" + postfix, -0.01, -1.5, 1.5);
  rooContainer->AddFormulaVar("CMS_hgg_modlin0" + postfix, "@0*@0", "CMS_hgg_lin0" + postfix);
  
  // prefix for models parameters
  std::map<int,std::string> parnames;
  parnames[1] = "modlin";
  parnames[2] = "modquad";
  parnames[3] = "modcubic";
  parnames[4] = "modquartic";
  parnames[5] = "modpol5_";
  parnames[6] = "modpol6_";
    
  // map order to categories flags + parameters names
  std::map<int, std::pair<std::vector<int>, std::vector<std::string> > > catmodels;
  // fill the map
  for(int icat=0; icat<nCategories_; ++icat) {
    // get the poly order for this category
    int catmodel = bkgPolOrderByCat[icat];
    std::vector<int> & catflags = catmodels[catmodel].first;
    std::vector<std::string> & catpars = catmodels[catmodel].second;
    // if this is the first time we find this order, build the parameters
    if( catflags.empty() ) {
      assert( catpars.empty() );
      // by default no category has the new model
      catflags.resize(nCategories_, 0);
      std::string & parname = parnames[catmodel];
      for(int iorder = 0; iorder<catmodel; ++iorder) {
	catpars.push_back( Form( "CMS_hgg_%s%d%s", parname.c_str(), iorder, +postfix.c_str() ) );
      }
    } else {
      assert( catflags.size() == nCategories_ && catpars.size() == catmodel );
    }
    // chose category order
    catflags[icat] = 1;
  }
  
  // now loop over the models and allocate the pdfs
  /// for(size_t imodel=0; imodel<catmodels.size(); ++imodel ) {
  for(std::map<int, std::pair<std::vector<int>, std::vector<std::string> > >::iterator modit = catmodels.begin();
      modit!=catmodels.end(); ++modit ) {
    std::vector<int> & catflags = modit->second.first;
    std::vector<std::string> & catpars = modit->second.second;
    
    rooContainer->AddSpecificCategoryPdf(&catflags[0],"data_pol_model" + postfix, "0", "CMS_hgg_mass", catpars, 70+catpars.size()); 
    // >= 71 means RooBernstein of order >= 1
  }
}

void initRooContainer(RooContainer* rooContainer, map<string, string> values, Smearing smearing) {

  rooContainer->SetNCategories(getUint(values, "nCategories"));
  rooContainer->nsigmas = getUint(values, "nSystSteps");
  rooContainer->sigmaRange = getUint(values, "systRange");
  
  smearing.effSmearPars.categoryType = "2CatR9_EBEE";
  smearing.effSmearPars.n_categories = 4;
  smearing.effSmearPars.efficiency_file = getString(values, "efficiencyFile");
  
  smearing.diPhoEffSmearPars.n_categories = 8;
  smearing.diPhoEffSmearPars.efficiency_file = getString(values, "efficiencyFile");
  
  if(getUint(values, "doEcorrectionSmear")) {
    smearing.photonSmearers_.push_back(smearing.eCorrSmearer);
  }

  if(getUint(values, "doEscaleSmear")) {
    if(getUint(values, "splitEscaleSyst")) {
      EnergySmearer::energySmearingParameters::eScaleVector tmp_scale_offset;
      EnergySmearer::energySmearingParameters::phoCatVector tmp_scale_cat;
      readEnergyScaleOffsets(getString(values, "scale_offset_corr_error_file"), tmp_scale_offset, tmp_scale_cat, false);
      assert(tmp_scale_offset.size() == 1); assert( ! tmp_scale_cat.empty() );
      
      smearing.eScaleCorrPars.categoryType = "Automagic";
      smearing.eScaleCorrPars.byRun = false;
      smearing.eScaleCorrPars.n_categories = tmp_scale_cat.size();
      smearing.eScaleCorrPars.photon_categories = tmp_scale_cat;
      
      smearing.eScaleCorrPars.scale_offset = tmp_scale_offset[0].scale_offset;
      smearing.eScaleCorrPars.scale_offset_error = tmp_scale_offset[0].scale_offset_error;
      
      smearing.eScaleCorrPars.smearing_sigma = tmp_scale_offset[0].scale_offset;
      smearing.eScaleCorrPars.smearing_sigma_error = tmp_scale_offset[0].scale_offset_error;
      
      EnergySmearer::energySmearingParameters::phoCatVectorIt icat = tmp_scale_cat.begin();
      for( ; icat != tmp_scale_cat.end(); ++icat ) {
	EnergySmearer * theSmear = new EnergySmearer(smearing.eScaleSmearer, EnergySmearer::energySmearingParameters::phoCatVector(1,*icat));
	theSmear->name(smearing.eScaleSmearer->name()+"_"+icat->name );
	theSmear->syst_only(true);
	std::cout << "Uncorrelated single photon category smearer " << theSmear->name() << std::endl; 
	smearing.photonSmearers_.push_back(theSmear);
	smearing.eScaleSmearers_.push_back(theSmear);
      }
      
      smearing.eScaleCorrSmearer = new EnergySmearer(smearing.eScaleCorrPars);
      smearing.eScaleCorrSmearer->name("E_scaleCorr");
      smearing.eScaleCorrSmearer->doEnergy(true);
      smearing.eScaleCorrSmearer->scaleOrSmear(true);
      smearing.eScaleCorrSmearer->syst_only(true);
      smearing.photonSmearers_.push_back(smearing.eScaleCorrSmearer);
      smearing.eScaleSmearers_.push_back(smearing.eScaleCorrSmearer);
      
      std::cout << "Uncorrelated single photon category smearer " << smearing.eScaleCorrSmearer->name() << std::endl; 
      
    } else {
      smearing.photonSmearers_.push_back(smearing.eScaleSmearer);
      smearing.eScaleSmearers_.push_back(smearing.eScaleSmearer);
    }
  }
  if(getUint(values, "doEresolSmear")) {
    if(getUint(values, "splitEresolSyst")) {

      // Use the same format used for the run-dependent energy corrections
      EnergySmearer::energySmearingParameters::eScaleVector tmp_smearing;
      EnergySmearer::energySmearingParameters::phoCatVector tmp_smearing_cat;
      readEnergyScaleOffsets(getString(values, "corr_smearing_file"), tmp_smearing, tmp_smearing_cat, false);
      
      // make sure that the scale correction and smearing info is as expected
      assert( tmp_smearing.size() == 1 );
      assert( ! tmp_smearing_cat.empty() );
      
      // copy the read info to the smarer parameters
      smearing.eResolCorrPars.categoryType = "Automagic";
      smearing.eResolCorrPars.byRun = false;
      smearing.eResolCorrPars.n_categories = tmp_smearing_cat.size();
      smearing.eResolCorrPars.photon_categories = tmp_smearing_cat;
      
      smearing.eResolCorrPars.scale_offset = tmp_smearing[0].scale_offset;
      smearing.eResolCorrPars.scale_offset_error = tmp_smearing[0].scale_offset_error;
      
      smearing.eResolCorrPars.smearing_sigma = tmp_smearing[0].scale_offset;
      smearing.eResolCorrPars.smearing_sigma_error = tmp_smearing[0].scale_offset_error;
      
      EnergySmearer::energySmearingParameters::phoCatVectorIt icat = tmp_smearing_cat.begin();
      for( ; icat != tmp_smearing_cat.end(); ++icat ) {
	EnergySmearer * theSmear = new EnergySmearer(smearing.eResolSmearer, EnergySmearer::energySmearingParameters::phoCatVector(1,*icat) );
	theSmear->name(smearing.eResolSmearer->name()+"_"+icat->name );
	std::cout << "Uncorrelated single photon category smearer " << theSmear->name() << std::endl; 
	smearing.photonSmearers_.push_back(theSmear);
	smearing.eResolSmearers_.push_back(theSmear);
      }
      
      smearing.eResolCorrSmearer = new EnergySmearer(smearing.eResolCorrPars);
      smearing.eResolCorrSmearer->name("E_resCorr");
      smearing.eResolCorrSmearer->doEnergy(true);
      smearing.eResolCorrSmearer->scaleOrSmear(true);
      smearing.eResolCorrSmearer->syst_only(true);
      smearing.photonSmearers_.push_back(smearing.eResolCorrSmearer);
      smearing.eResolSmearers_.push_back(smearing.eResolCorrSmearer);
      std::cout << "Uncorrelated single photon category smearer " << smearing.eResolCorrSmearer->name() << std::endl; 
      
    } else {
      smearing.photonSmearers_.push_back(smearing.eResolSmearer);
      smearing.eResolSmearers_.push_back(smearing.eResolSmearer);
    }
  }
  if(getUint(values, "doPhotonIdEffSmear")) {
    // photon ID efficiency 
    std::cerr << __LINE__ << std::endl; 
    smearing.idEffSmearer = new EfficiencySmearer(smearing.effSmearPars);
    smearing.idEffSmearer->name("idEff");
    smearing.idEffSmearer->setEffName("ratioTP");
    smearing.idEffSmearer->init();
    smearing.idEffSmearer->doPhoId(true);
    smearing.photonSmearers_.push_back(smearing.idEffSmearer);
  }
  if(getUint(values, "doR9Smear")) {
    // R9 re-weighting
    smearing.r9Smearer = new EfficiencySmearer(smearing.effSmearPars);
    smearing.r9Smearer->name("r9Eff");
    smearing.r9Smearer->setEffName("ratioR9");
    smearing.r9Smearer->init();
    smearing.r9Smearer->doR9(true);
    smearing.photonSmearers_.push_back(smearing.r9Smearer);
  }
  if(getUint(values, "doVtxEffSmear")) {
    // Vertex ID
    std::cerr << __LINE__ << std::endl; 
    smearing.vtxEffSmearer = new DiPhoEfficiencySmearer(smearing.diPhoEffSmearPars);   // triplicate TF1's here
    smearing.vtxEffSmearer->name("vtxEff");
    smearing.vtxEffSmearer->setEffName("ratioVertex");
    smearing.vtxEffSmearer->doVtxEff(true);
    smearing.vtxEffSmearer->init();
    smearing.diPhotonSmearers_.push_back(smearing.vtxEffSmearer);
  }
  if(getUint(values, "doTriggerEffSmear")) {
    // trigger efficiency
    std::cerr << __LINE__ << std::endl; 
    smearing.triggerEffSmearer = new DiPhoEfficiencySmearer(smearing.diPhoEffSmearPars);
    smearing.triggerEffSmearer->name("triggerEff");
    smearing.triggerEffSmearer->setEffName("effL1HLT");
    smearing.triggerEffSmearer->doVtxEff(false);
    smearing.triggerEffSmearer->init();
    smearing.diPhotonSmearers_.push_back(smearing.triggerEffSmearer);
  }
  if(getUint(values, "doKFactorSmear")) {
    // kFactor efficiency
    std::cerr << __LINE__ << std::endl; 
    smearing.kFactorSmearer = new KFactorSmearer(getString(values, "kfacHist"));
    smearing.kFactorSmearer->name("kFactor");
    smearing.kFactorSmearer->init();
    smearing.genLevelSmearers_.push_back(smearing.kFactorSmearer);
  }
  if(getUint(values, "doInterferenceSmear")) {
    // interference efficiency
    std::cerr << __LINE__ << std::endl; 
    smearing.interferenceSmearer = new InterferenceSmearer(2.5e-2, 0.);
    smearing.genLevelSmearers_.push_back(smearing.interferenceSmearer);
  }
  
  // Define the number of categories for the statistical analysis and
  // the systematic sets to be formed
  
  if(getUint(values, "doEcorrectionSmear") && getUint(values, "doEcorrectionSyst")) {
    // instance of this smearer done in PhotonAnalysis
    smearing.systPhotonSmearers_.push_back(smearing.eCorrSmearer);
    std::vector<std::string> sys(1, smearing.eCorrSmearer->name());
    std::vector<int> sys_t(1, -1);   // -1 for signal, 1 for background 0 for both
    rooContainer->MakeSystematicStudy(sys,sys_t);
  }

  if(getUint(values, "doEscaleSmear") && getUint(values, "doEscaleSyst")) {
    std::cout << smearing.eScaleSmearers_.size() << std::endl;
    for(std::vector<EnergySmearer*>::iterator ei=smearing.eScaleSmearers_.begin(); ei!=smearing.eScaleSmearers_.end(); ++ei){
      smearing.systPhotonSmearers_.push_back( *ei );
      std::vector<std::string> sys(1,(*ei)->name());
      std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
      rooContainer->MakeSystematicStudy(sys,sys_t);
    }
  }

  if(getUint(values, "doEresolSmear") && getUint(values, "doEresolSyst")) {
    for(std::vector<EnergySmearer*>::iterator ei=smearing.eResolSmearers_.begin(); ei!=smearing.eResolSmearers_.end(); ++ei){
      smearing.systPhotonSmearers_.push_back( *ei );
      std::vector<std::string> sys(1,(*ei)->name());
      std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
      rooContainer->MakeSystematicStudy(sys,sys_t);
    }
  }

  if(getUint(values, "doPhotonIdEffSmear") && getUint(values, "doPhotonIdEffSyst")) {
    smearing.systPhotonSmearers_.push_back(smearing.idEffSmearer);
    std::vector<std::string> sys(1, smearing.idEffSmearer->name());
    std::vector<int> sys_t(1, -1);   // -1 for signal, 1 for background 0 for both
    rooContainer->MakeSystematicStudy(sys, sys_t);
  }
  if(getUint(values, "doR9Smear") && getUint(values, "doR9Syst")) {
    smearing.systPhotonSmearers_.push_back(  smearing.r9Smearer);
    std::vector<std::string> sys(1,smearing.r9Smearer->name());
    std::vector<int> sys_t(1, -1);   // -1 for signal, 1 for background 0 for both
    rooContainer->MakeSystematicStudy(sys, sys_t);
  }
  if(getUint(values, "doVtxEffSmear") && getUint(values, "doVtxEffSyst")) {
    smearing.systDiPhotonSmearers_.push_back(smearing.vtxEffSmearer);
    std::vector<std::string> sys(1, smearing.vtxEffSmearer->name());
    std::vector<int> sys_t(1, -1);   // -1 for signal, 1 for background 0 for both
    rooContainer->MakeSystematicStudy(sys, sys_t);
  }
  if(getUint(values, "doTriggerEffSmear") && getUint(values, "doTriggerEffSyst")) {
    smearing.systDiPhotonSmearers_.push_back(smearing.triggerEffSmearer);
    std::vector<std::string> sys(1, smearing.triggerEffSmearer->name());
    std::vector<int> sys_t(1, -1);   // -1 for signal, 1 for background 0 for both
    rooContainer->MakeSystematicStudy(sys, sys_t);
  }
  if(getUint(values, "doKFactorSmear") && getUint(values, "doKFactorSyst")) {
    smearing.systGenLevelSmearers_.push_back(smearing.kFactorSmearer);
    std::vector<std::string> sys(1, smearing.kFactorSmearer->name());
    std::vector<int> sys_t(1, -1);   // -1 for signal, 1 for background 0 for both
    rooContainer->MakeSystematicStudy(sys, sys_t);
  }
  
  rooContainer->AddGlobalSystematic("lumi", 1.045, 1.00);
  
  // Create observables for shape-analysis with ranges
  rooContainer->AddObservable("CMS_hgg_mass", getFloat(values, "massMin"), getFloat(values, "massMax"));
  rooContainer->AddConstant("IntLumi", getFloat(values, "intlumi"));

  // SM Model
  rooContainer->AddConstant("XSBR_tth_155", 0.00004370);
  rooContainer->AddConstant("XSBR_ggh_150", 0.01428);
  rooContainer->AddConstant("XSBR_vbf_150", 0.001308);
  rooContainer->AddConstant("XSBR_wzh_150", 0.000641);
  rooContainer->AddConstant("XSBR_tth_150", 0.000066);
  rooContainer->AddConstant("XSBR_ggh_145", 0.018820);
  rooContainer->AddConstant("XSBR_vbf_145", 0.001676);
  rooContainer->AddConstant("XSBR_wzh_145", 0.000891);
  rooContainer->AddConstant("XSBR_tth_145", 0.000090);
  rooContainer->AddConstant("XSBR_ggh_140", 0.0234109);
  rooContainer->AddConstant("XSBR_vbf_140", 0.00203036);
  rooContainer->AddConstant("XSBR_wzh_140", 0.001163597);
  rooContainer->AddConstant("XSBR_tth_140", 0.000117189);
  rooContainer->AddConstant("XSBR_ggh_135", 0.0278604);
  rooContainer->AddConstant("XSBR_vbf_135", 0.002343);
  rooContainer->AddConstant("XSBR_wzh_135", 0.001457559);
  rooContainer->AddConstant("XSBR_tth_135", 0.000145053);
  rooContainer->AddConstant("XSBR_ggh_130", 0.0319112);
  rooContainer->AddConstant("XSBR_vbf_130", 0.00260804);
  rooContainer->AddConstant("XSBR_wzh_130", 0.001759636);
  rooContainer->AddConstant("XSBR_tth_130", 0.000173070);
  rooContainer->AddConstant("XSBR_ggh_125", 0.0350599);
  rooContainer->AddConstant("XSBR_vbf_125", 0.00277319);
  rooContainer->AddConstant("XSBR_wzh_125", 0.002035123);
  rooContainer->AddConstant("XSBR_tth_125", 0.000197718);
  rooContainer->AddConstant("XSBR_ggh_120", 0.0374175);
  rooContainer->AddConstant("XSBR_vbf_120", 0.00285525);
  rooContainer->AddConstant("XSBR_wzh_120", 0.002285775);
  rooContainer->AddConstant("XSBR_tth_120", 0.00021951);
  rooContainer->AddConstant("XSBR_ggh_123", 0.0360696);
  rooContainer->AddConstant("XSBR_vbf_123", 0.00281352);
  rooContainer->AddConstant("XSBR_wzh_123", 0.00213681);
  rooContainer->AddConstant("XSBR_tth_123", 0.00020663);
  rooContainer->AddConstant("XSBR_ggh_121", 0.0369736);
  rooContainer->AddConstant("XSBR_vbf_121", 0.00284082);
  rooContainer->AddConstant("XSBR_wzh_121", 0.00223491);
  rooContainer->AddConstant("XSBR_tth_121", 0.00021510);
  rooContainer->AddConstant("XSBR_ggh_115", 0.0386169);
  rooContainer->AddConstant("XSBR_vbf_115", 0.00283716);
  rooContainer->AddConstant("XSBR_wzh_115", 0.002482089);
  rooContainer->AddConstant("XSBR_tth_115", 0.000235578);
  rooContainer->AddConstant("XSBR_ggh_110", 0.0390848);
  rooContainer->AddConstant("XSBR_vbf_110", 0.00275406);
  rooContainer->AddConstant("XSBR_wzh_110", 0.002654575);
  rooContainer->AddConstant("XSBR_tth_110", 0.000247629);
  rooContainer->AddConstant("XSBR_ggh_105", 0.0387684);
  rooContainer->AddConstant("XSBR_vbf_105", 0.00262016);
  rooContainer->AddConstant("XSBR_wzh_105", 0.002781962);
  rooContainer->AddConstant("XSBR_tth_105", 0.000255074);
  
  // FF model 
  rooContainer->AddConstant("ff_XSBR_vbf_150", 0.00259659);
  rooContainer->AddConstant("ff_XSBR_wzh_150", 0.00127278);
  rooContainer->AddConstant("ff_XSBR_vbf_145", 0.00387544);
  rooContainer->AddConstant("ff_XSBR_wzh_145", 0.00205969);
  rooContainer->AddConstant("ff_XSBR_vbf_140", 0.00565976);
  rooContainer->AddConstant("ff_XSBR_wzh_140", 0.003243602);
  rooContainer->AddConstant("ff_XSBR_vbf_135", 0.00825);
  rooContainer->AddConstant("ff_XSBR_wzh_135", 0.00513225);
  rooContainer->AddConstant("ff_XSBR_vbf_130", 0.0122324);
  rooContainer->AddConstant("ff_XSBR_wzh_130", 0.00825316);
  rooContainer->AddConstant("ff_XSBR_vbf_125", 0.0186494);
  rooContainer->AddConstant("ff_XSBR_wzh_125", 0.01368598);
  rooContainer->AddConstant("ff_XSBR_vbf_123", 0.022212);
  rooContainer->AddConstant("ff_XSBR_wzh_123", 0.0168696);
  rooContainer->AddConstant("ff_XSBR_vbf_121", 0.0266484);
  rooContainer->AddConstant("ff_XSBR_wzh_121", 0.0209646);
  rooContainer->AddConstant("ff_XSBR_vbf_120", 0.0293139);
  rooContainer->AddConstant("ff_XSBR_wzh_120", 0.02346729);
  rooContainer->AddConstant("ff_XSBR_vbf_115", 0.0482184);
  rooContainer->AddConstant("ff_XSBR_wzh_115", 0.04218386);
  rooContainer->AddConstant("ff_XSBR_vbf_110", 0.083181);
  rooContainer->AddConstant("ff_XSBR_wzh_110", 0.08017625);
  rooContainer->AddConstant("ff_XSBR_vbf_105", 0.151616);
  rooContainer->AddConstant("ff_XSBR_wzh_105", 0.1609787);
  
  // -----------------------------------------------------
  // Configurable background model
  // if no configuration was given, set some defaults

  int nCategories_ = getUint(values, "nCategories");
  int nInclusiveCategories_ = getUint(values, "nInclusiveCategories");
  int nVBFCategories_ = getUint(values, "nVBFCategories");
  int nVHhadCategories_ = getUint(values, "nVHhadCategories");
  int nVHlepCategories_ = getUint(values, "nVHlepCategories");
  std::vector<int> bkgPolOrderByCat = getVint(values, "bkgPolOrderByCat");

  if( bkgPolOrderByCat.empty() ) {
    for(int i=0; i<nCategories_; i++){
      if(i<nInclusiveCategories_) {
	bkgPolOrderByCat.push_back(5);
      } else if(i<nInclusiveCategories_+nVBFCategories_){
	bkgPolOrderByCat.push_back(3);
      } else if(i<nInclusiveCategories_+nVBFCategories_+nVHhadCategories_){
	bkgPolOrderByCat.push_back(2);
      } else if(i<nInclusiveCategories_+nVBFCategories_+nVHhadCategories_+nVHlepCategories_){
	bkgPolOrderByCat.push_back(1);
      }
    }
  }

  // build the model
  buildBkgModel(rooContainer, values);
  
  /// -----------------------------------
  // -----------------------------------------------------
  // Make some data sets from the observables to fill in the event loop         
  // Binning is for histograms (will also produce unbinned data sets)
  unsigned int nDataBins = getUint(values, "nDataBins");
  rooContainer->CreateDataSet("CMS_hgg_mass", "data_mass", nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
  rooContainer->CreateDataSet("CMS_hgg_mass", "bkg_mass" , nDataBins);            
  

  std::vector<int> sigPointsToBook = getVint(values, "signalPoints");
  // Create Signal DataSets:
  for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
    int sig = sigPointsToBook[isig];
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_ggh_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_vbf_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_wzh_mass_m%d", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_tth_mass_m%d", sig), nDataBins);   
    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_ggh_mass_m%d_rv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_vbf_mass_m%d_rv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_wzh_mass_m%d_rv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_tth_mass_m%d_rv", sig), nDataBins);    
    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_ggh_mass_m%d_wv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_vbf_mass_m%d_wv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_wzh_mass_m%d_wv", sig), nDataBins);    
    rooContainer->CreateDataSet("CMS_hgg_mass", Form("sig_tth_mass_m%d_wv", sig), nDataBins);    
  }
  
  // Make more datasets representing Systematic Shifts of various quantities
  for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
    int sig = sigPointsToBook[isig];
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_ggh_mass_m%d", sig), -1);    
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_vbf_mass_m%d", sig), -1);    
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_wzh_mass_m%d", sig), -1);    
    rooContainer->MakeSystematics("CMS_hgg_mass", Form("sig_tth_mass_m%d", sig), -1);    
  }
}

void fitToData(RooContainer* rooContainer, std::map<std::string, std::string> values) {

  std::string postfix = (getUint(values,"dataIs2011")?"":"_8TeV");
  rooContainer->FitToData("data_pol_model" + postfix, "data_mass");  // Fit to full range of dataset
}

//----------------------------------------------------------------------
void initConfig(std::map<std::string, std::string> values, HistoContainer *histoContainer, RooContainer* rooContainer, Smearing smearing) {

  // FIXME ADD A CONFIGURABLE SWITCH;
  initRooContainer(rooContainer, values, smearing); 
  bookHistogram(histoContainer, values);
}

//----------------------------------------------------------------------


GenericAnalysis *openAnalysisCode(const string &fname)
{
  // open the shared library
  void* soHandle = dlopen(fname.c_str(), RTLD_LAZY);

  if (! soHandle)
    {
      cerr << "could not open analysis shared object file '" << fname << "': " << dlerror() << endl;
      exit(1);
    }

    // load the instantiation function (the function which creates
    // an instance of the needed class)
    typedef GenericAnalysis *(*InstantiationFunctionType)();

    dlerror();

    const char *instantiationFunctionName = "makeInstance";

    InstantiationFunctionType instantiationFunction = (InstantiationFunctionType) dlsym(soHandle, instantiationFunctionName);
    const char *errMessage = dlerror();
    if (errMessage)
    {
      cerr << "could not find symbol '" << instantiationFunctionName << "' in analysis shared object file '" << fname << "'. Exiting." << endl;
      dlclose(soHandle);
      exit(1);
    }

    GenericAnalysis *analysis = instantiationFunction();
    return analysis;

}


//----------------------------------------------------------------------

void usage()
{
  cerr << endl
       << "usage:   plotter plotvariables.dat myanalysis.so input.root output.root" << endl
       << endl
    ;
  exit(1);
}

//----------------------------------------------------------------------
int main(int argc, char **argv)
{
  if (argc != 4 + 1)
    usage();

  string configFname = argv[1];
  string analysisCodeFname = argv[2];
  string inputFname = argv[3];
  string outputFname = argv[4];

  HistoContainer *histoContainer = new HistoContainer();
  RooContainer *rooContainer = new RooContainer();
  Smearing smearing;

  // read the configuration file and book histograms 
  map<std::string, std::string> values = parseConfigFile(configFname);
  initConfig(values, histoContainer, rooContainer, smearing);

  //--------------------

  TRint *myapp = NULL;
  if (interactive)
  {
    int dummyArgc = 1;
    char* dummyArgv[2];
    dummyArgv[0] = argv[0];
    dummyArgv[1] = NULL;

    myapp = new TRint("rint", &dummyArgc, &dummyArgv[0]);
  }

  //--------------------
  // open the input file
  //--------------------

  TFile *fin = TFile::Open(inputFname.c_str());

  if (fin == NULL || !fin->IsOpen())
  {
    cerr << "could not open input file '" << inputFname << "'" << endl;
    exit(1);
  }

  TTree *tree = (TTree*)(fin->Get(inputTreeName.c_str()));
  assert(tree != NULL);

  //--------------------
  // open the analysis shared object
  //--------------------
  GenericAnalysis *analysis = openAnalysisCode(analysisCodeFname);

  // FIXME HOW TO PASS PARAMETER TO THE ANALYSIS ???
  //analysis.setParamters();

  //--------------------

  // must also activate branches
  analysis->setBranchAddresses(tree, values);

  // loop on the events
  unsigned numEvents = tree->GetEntries();
  for (unsigned i = 0; i < numEvents; ++i)
  {
    // read ith event
    tree->GetEntry(i);

    // call user function
    analysis->analyze(histoContainer, rooContainer, smearing);

  } // end of loop over events

  //--------------------
  
  fitToData(rooContainer, values);

  TFile *fout = TFile::Open(outputFname.c_str(), "RECREATE");

  // write out histograms to output file
  fout->cd();
  histoContainer->Save();
  rooContainer->Save();

  // do plotting FIXME
  Plotting plotter(histoContainer);

  // see http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/Macros/Normalization_8TeV.cc?revision=1.10&view=markup
  // for itype numbers for each mass

  // mh = 125 GeV
  plotter.addSignalItype(-37);
  plotter.addSignalItype(-38);
  plotter.addSignalItype(-39);
  plotter.addSignalItype(-40);
  
  plotter.plotAll();

  if (interactive)
  {
    myapp->Run(kTRUE);            // run ROOT interactively
  }

  // close output file
  fout->Close();

}
