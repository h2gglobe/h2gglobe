#include "../interface/PhotonAnalysis.h"


#include "PhotonReducedInfo.h"
#include "Sorters.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#define PADEBUG 0

using namespace std;



// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::PhotonAnalysis()  : 
	runStatAnalysis(false), doTriggerSelection(false),
	name_("PhotonAnalysis"),
	vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams),
	tmvaPerVtxMethod("BDTG"),
	tmvaPerVtxWeights(""),
	tmvaPerEvtMethod("evtBTG"),
	tmvaPerEvtWeights(""),
	energyCorrectionMethod("DaunceyAndKenzie"), energyCorrected(0), energyCorrectedError(0)
{
	addConversionToMva=true;
	mvaVertexSelection=false;
	useDefaultVertex=false;
	forcedRho = -1.;

	keepPP = true;
	keepPF = true;
	keepFF = true; 
}

// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::~PhotonAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void readEnergyScaleOffsets(const std::string &fname, EnergySmearer::energySmearingParameters::eScaleVector &escaleOffsets)
{
  // read in energy scale corrections to be applied in run ranges
	std::fstream in(fname.c_str());
	assert( in );
	char line[200];
	float EBHighR9, EBLowR9, EBm4HighR9, EBm4LowR9, EEHighR9, EELowR9; 
	int  first, last;
	do {
		in.getline( line, 200, '\n' );

		if( sscanf(line,"%d %d %f %f %f %f",&first, &last, &EBHighR9, &EBLowR9, &EBm4HighR9, &EBm4LowR9, &EEHighR9, &EELowR9) != 6 ) { continue; } 
		std::cerr << "Energy scale by run " <<  first<< " " <<  last<< " " <<  EBHighR9<< " " <<  EBLowR9 << " " <<  EBm4HighR9<< " " <<  EBm4LowR9<< " " <<  EEHighR9<< " " <<  EELowR9 << std::endl;

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

	} while( in );
	
	in.close();
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuMap(const char * fname, TDirectory * dir)
{
	std::fstream in(fname);
	assert( in );
	char line[200];
	int typid;
	char dname[200];
	do {
		in.getline( line, 200, '\n' );

		if( sscanf(line,"%d %s",&typid,dname) != 2 ) { continue; } 
		std::cerr << "Reading PU weights for sample " << typid << " from " << dname << std::endl;
		TDirectory * subdir = (TDirectory *)dir->Get(dname);
		assert( subdir != 0 );
		loadPuWeights(typid, subdir);
	} while ( in );
	in.close();
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuWeights(int typid, TDirectory * puFile)
{
            cout<<"Reweighting events for pileup."<<endl;
	    TH1 * hweigh = (TH1*) puFile->Get("weights");
	    if( hweigh == 0 ) {
		    hweigh = (TH1*) puFile->Get("NPUWeights");
	    }
	    if( hweigh != 0 ) { 
		    cout<< " This is a pre-processed pileup reweighing file." <<endl;
		    TH1 * gen_pu = (TH1*)puFile->Get("generated_pu");
		    if( gen_pu == 0 ) {
			    gen_pu = (TH1*)puFile->Get("NPUSource");
		    }
		    // Normalize weights such that the total cross section is unchanged
		    TH1 * eff = (TH1*)hweigh->Clone("eff");
		    eff->Multiply(gen_pu);
		    hweigh->Scale( gen_pu->Integral() / eff->Integral()  );
		    weights[typid].clear();
		    for( int ii=1; ii<hweigh->GetNbinsX(); ++ii ) {
			    weights[typid].push_back(hweigh->GetBinContent(ii)); 
		    }
	    } 
	    //// else {
	    //// 	    TH1D * histo = (TH1D*) puFile->Get("pileup");
	    //// 	    if( histo != 0 ) {
	    //// 		    weights[typid] = l.generate_flat10_weights(histo);
	    //// 	    }
	    //// }
	    std::cout << "pile-up weights: ["<<typid<<"]";
	    std::copy(weights[typid].begin(), weights[typid].end(), std::ostream_iterator<double>(std::cout,","));
	    std::cout << std::endl;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Init(LoopAll& l) 
{
	if(PADEBUG) 
		cout << "InitRealPhotonAnalysis START"<<endl;

	if(energyCorrectionMethod=="DaunceyAndKenzie"){
		energyCorrected		= (l.pho_residCorrEnergy);
		energyCorrectedError= (l.pho_residCorrResn);
	}else if(energyCorrectionMethod=="Bendavid"){
		energyCorrected		= (l.pho_regr_energy);
		energyCorrectedError= (l.pho_regr_energyerr);
//	}else if(energyCorrectionMethod=="PFRegression"){
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
	triggerSelections.push_back(TriggerSelection(161216,165087));
	triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
	triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
	triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
	triggerSelections.back().addpath("HLT_Photon20_R9Id_Photon18_R9Id_v");

        // /cdaq/physics/Run2011/1e33/v1.3/HLT/V2
	triggerSelections.push_back(TriggerSelection(165088,173198));
	triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
	triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
	triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
	triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");

        // /cdaq/physics/Run2011/3e33/v1.1/HLT/V1
	triggerSelections.push_back(TriggerSelection(173236,178380));
	triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v");
	triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v");
	triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v");
	triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
	
        // /cdaq/physics/Run2011/5e33/v1.4/HLT/V3
	triggerSelections.push_back(TriggerSelection(178420,-1));
	triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v");
	triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v");
	triggerSelections.back().addpath("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v");
	triggerSelections.back().addpath("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v");
	
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
		l.SetPhotonCutsInCategories((LoopAll::phoCiCIDLevel)iLevel, &cic6_cuts_lead[0][0], &cic6_cuts_sublead[0][0], &cic4_cuts_lead[0][0], &cic4_cuts_sublead[0][0] );
		
		float * cic6_cuts_arrays_lead[phoNCUTS] = {
			&l.cic6_cut_lead_isosumoet[0][0], &l.cic6_cut_lead_isosumoetbad[0][0], &l.cic6_cut_lead_trkisooet[0][0], &l.cic6_cut_lead_sieie[0][0],
			&l.cic6_cut_lead_hovere[0][0], &l.cic6_cut_lead_r9[0][0], &l.cic6_cut_lead_drtotk_25_99[0][0], &l.cic6_cut_lead_pixel[0][0] 
		};
		
		float * cic6_cuts_arrays_sublead[phoNCUTS] = {
			&l.cic6_cut_sublead_isosumoet[0][0], &l.cic6_cut_sublead_isosumoetbad[0][0], &l.cic6_cut_sublead_trkisooet[0][0], 
			&l.cic6_cut_sublead_sieie[0][0], &l.cic6_cut_sublead_hovere[0][0], &l.cic6_cut_sublead_r9[0][0],
			&l.cic6_cut_sublead_drtotk_25_99[0][0], &l.cic6_cut_sublead_pixel[0][0]
		};

		float * cic4_cuts_arrays_lead[phoNCUTS] = {
			&l.cic4_cut_lead_isosumoet[0][0], &l.cic4_cut_lead_isosumoetbad[0][0], &l.cic4_cut_lead_trkisooet[0][0], &l.cic4_cut_lead_sieie[0][0],
			&l.cic4_cut_lead_hovere[0][0], &l.cic4_cut_lead_r9[0][0], &l.cic4_cut_lead_drtotk_25_99[0][0], &l.cic4_cut_lead_pixel[0][0] 
		};
		
		float * cic4_cuts_arrays_sublead[phoNCUTS] = {
			&l.cic4_cut_sublead_isosumoet[0][0], &l.cic4_cut_sublead_isosumoetbad[0][0], &l.cic4_cut_sublead_trkisooet[0][0], 
			&l.cic4_cut_sublead_sieie[0][0], &l.cic4_cut_sublead_hovere[0][0], &l.cic4_cut_sublead_r9[0][0],
			&l.cic4_cut_sublead_drtotk_25_99[0][0], &l.cic4_cut_sublead_pixel[0][0]
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
		}
	}
	
	if( tmvaPerVtxWeights != ""  ) {
		tmvaPerVtxVariables_.push_back("ptbal"), tmvaPerVtxVariables_.push_back("ptasym"), tmvaPerVtxVariables_.push_back("logsumpt2");
		if( addConversionToMva ) {
			tmvaPerVtxVariables_.push_back("limPullToConv");
			tmvaPerVtxVariables_.push_back("nConv");
		}
		tmvaPerVtxReader_ = new TMVA::Reader( "!Color:!Silent" );
		HggVertexAnalyzer::bookVariables( *tmvaPerVtxReader_, tmvaPerVtxVariables_ );
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
	eSmearDataPars.n_categories = 6;
	std::cerr << "Reading energy scale offsets " << scale_offset_file << std::endl;
	readEnergyScaleOffsets(scale_offset_file, eSmearDataPars.scale_offset_byrun);
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
	
	// energy scale systematics to MC
        eScaleSmearer = new EnergySmearer( eSmearPars );
	eScaleSmearer->name("E_scale");
	eScaleSmearer->doEnergy(true);
	eScaleSmearer->scaleOrSmear(true);

	if( doEcorrectionSmear ) {
	  eCorrSmearer = new EnergySmearer( eSmearPars );
	  eCorrSmearer->name("E_corr");
	  // activating pho corrections to this instance of EnergySmearer, implies that it won't touch Escale and Eresolution
	  eCorrSmearer->doCorrections(true); 
	}
	
	if (l.typerun == 2 || l.typerun == 1) {
	}
	
    /* -------------------------------------------------------------------------------------------
    Pileup Reweighting
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
    ----------------------------------------------------------------------------------------------  */
    if (puHist != "") {
        if(PADEBUG) 
		cout << "Opening PU file"<<endl;
        TFile* puFile = TFile::Open( puHist );
        if (puFile) {
		if( puMap != "" ) {
		  loadPuMap(puMap, puFile); 
		} else {
		  loadPuWeights(0, puFile);
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
    } 

	if(PADEBUG) 
		cout << "InitRealPhotonAnalysis END"<<endl;

	// FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
        if(PADEBUG) 
                cout << "Analysis START"<<endl;
        pho_presel.clear();

        //remove process ID 18 from gamma+jet to avoid double counting with born+box
        if (l.itype[l.current]==3 && l.process_id==18) return;

	//apply pileup reweighting
	unsigned int n_pu = l.pu_n;
	float weight =1.;
	if (l.itype[l.current] !=0 && puHist != "") {
	  std::vector<double> & puweights = weights.find( l.itype[l.current] ) != weights.end() ? weights[ l.itype[l.current] ] : weights[0]; 
	  if(n_pu<puweights.size()){
	    weight *= puweights[n_pu];
	  }    
	  else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	    cout<<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<"), event will not be reweighted for pileup"<<endl;
	  }
	}

        if( pho_presel.empty() ) { 
                PreselectPhotons(l,jentry);
        }

	if( pho_acc.size() < 2 || pho_et[ pho_acc[0] ] < presel_scet1 ) return;

        int leadLevel=LoopAll::phoSUPERTIGHT, subLevel=LoopAll::phoSUPERTIGHT;

	//Fill histograms to use as denominator (pre-selection only) and numerator (selection applied)
	//for photon ID efficiency calculation.  To avoid ambiguities concerning vertex choice, use only 
	//events with one diphoton pair (close to 100% of signal events)
        if (l.dipho_n==1) {

	  int ivtx = l.dipho_vtxind[0];
	  int lead = l.dipho_leadind[0];
	  int sublead = l.dipho_subleadind[0];

	  TLorentzVector lead_p4 = l.get_pho_p4(lead,ivtx,&corrected_pho_energy[0]); 
	  TLorentzVector sublead_p4 = l.get_pho_p4(sublead,ivtx,&corrected_pho_energy[0]); 
	  float leadEta = ((TVector3 *)l.sc_xyz->At(l.pho_scind[lead]))->Eta();
	  float subleadEta = ((TVector3 *)l.sc_xyz->At(l.pho_scind[sublead]))->Eta();

	  //apply pre-selection
	  bool passpresel = true;
	  if(lead_p4.Pt() < 40. || sublead_p4.Pt() < 30. || 
	      fabs(leadEta) > 2.5 || fabs(subleadEta) > 2.5 || 
	      ( fabs(leadEta) > 1.4442 && fabs(leadEta) < 1.566 ) || ( fabs(subleadEta) > 1.4442 && fabs(subleadEta) < 1.566 ))
	    passpresel = false;
	  if (lead != sublead && passpresel) {

	    int leadpho_category = l.PhotonCategory(lead, 2, 2);
	    int subleadpho_category = l.PhotonCategory(sublead, 2, 2);

	    //Fill eta and pt distributions after pre-selection only (efficiency denominator)
	    l.FillHist("pho1_pt_presel",0,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_presel",0,sublead_p4.Pt(), weight);
	    l.FillHist("pho1_eta_presel",0,leadEta, weight);
	    l.FillHist("pho2_eta_presel",0,subleadEta, weight);

	    l.FillHist("pho1_pt_presel",leadpho_category+1,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_presel",subleadpho_category+1,sublead_p4.Pt(), weight);
	    l.FillHist("pho1_eta_presel",leadpho_category+1,leadEta, weight);
	    l.FillHist("pho2_eta_presel",subleadpho_category+1,subleadEta, weight);

	    //Apply single photon CiC selection and fill eta and pt distributions (efficiency numerator)
	    std::vector<std::vector<bool> > ph_passcut;
	    if( l.PhotonCiCSelectionLevel(lead, ivtx, ph_passcut, 4, 0, &corrected_pho_energy[0]) >=  (LoopAll::phoCiCIDLevel) leadLevel) {
	      l.FillHist("pho1_pt_sel",0,lead_p4.Pt(), weight);
	      l.FillHist("pho1_eta_sel",0,leadEta, weight);
	      l.FillHist("pho1_pt_sel",leadpho_category+1,lead_p4.Pt(), weight);
	      l.FillHist("pho1_eta_sel",leadpho_category+1,leadEta, weight);
	    }
	    if( l.PhotonCiCSelectionLevel(sublead, ivtx, ph_passcut, 4, 1, &corrected_pho_energy[0]) >=  (LoopAll::phoCiCIDLevel) subLevel ) {
	      l.FillHist("pho2_pt_sel",0,sublead_p4.Pt(), weight);
	      l.FillHist("pho2_eta_sel",0,subleadEta, weight);
	      l.FillHist("pho2_pt_sel",subleadpho_category+1,sublead_p4.Pt(), weight);
	      l.FillHist("pho2_eta_sel",subleadpho_category+1,subleadEta, weight);
	    }
	  }
        }

	//Apply diphoton CiC selection
	int dipho_id = l.DiphotonCiCSelection((LoopAll::phoCiCIDLevel) leadLevel, (LoopAll::phoCiCIDLevel) subLevel, 40., 30.0, 4, applyPtoverM, &corrected_pho_energy[0]);

        if (dipho_id > -1){
	  std::pair<int,int> diphoton_index = std::make_pair(l.dipho_leadind[dipho_id],l.dipho_subleadind[dipho_id]);
	  TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[dipho_id], l.dipho_vtxind[dipho_id], &corrected_pho_energy[0]);
	  TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[dipho_id], l.dipho_vtxind[dipho_id], &corrected_pho_energy[0]);
	  TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	  
          //// TLorentzVector *lead_p4 = (TLorentzVector*)l.pho_p4->At(diphoton_index.first);
          //// TLorentzVector *sublead_p4 = (TLorentzVector*)l.pho_p4->At(diphoton_index.second);
          //// TLorentzVector Higgs = *lead_p4 + *sublead_p4;

	  // Calculate the Rest Frame Quantities:
	  TLorentzVector *rest_lead_p4 = (TLorentzVector*) lead_p4.Clone();
	  TLorentzVector *rest_sublead_p4 = (TLorentzVector*) sublead_p4.Clone();
	  TVector3 higgsBoostVector = Higgs.BoostVector();
	  rest_lead_p4->Boost(higgsBoostVector);        // lead photon in Higgs Decay Frame
	  rest_sublead_p4->Boost(higgsBoostVector);    // sub-lead photon in Higgs Decay Frame

	  // Calculate variables to be filled
	  float decayAngle  = l.DeltaPhi(lead_p4.Phi(),sublead_p4.Phi());
	  float helicityAngle = fabs(lead_p4.E()-sublead_p4.E())/Higgs.P();
	  float maxeta = lead_p4.Eta();
	  if (fabs(sublead_p4.Eta())>fabs(lead_p4.Eta())) maxeta = sublead_p4.Eta();
	  float HT = lead_p4.Pt() + sublead_p4.Pt();
	  float mH = Higgs.M();

	  /*
	  //Alternative definition of background categories for QCD,GJET, requiring matching of gen
	  //photons to reco photons (main background category definition, based on gen info only is defined in
	  //PhotonAnalysis::SkimEvents)
	  BkgCategory bkgCat = promptprompt;
          if (l.itype[l.current]==1 || l.itype[l.current]==2 || l.itype[l.current]==3) {

            bool lead_is_prompt=false;
            bool sublead_is_prompt=false;

            for(int igp=0; igp<l.gp_n; ++igp) {

              int id = l.gp_pdgid[igp];
              int mother = l.gp_mother[igp];
              int motherid = l.gp_pdgid[mother];

              if ( id==22 && mother>-1 && (motherid==21 || abs(motherid)<=8 ) ) {
                //cout << "mother " << mother << " " << motherid << endl;
                TLorentzVector *gp_p4 = (TLorentzVector*)l.gp_p4->At(igp);
                if ( gp_p4->DeltaR(lead_p4) < 0.3) lead_is_prompt = true;
                if ( gp_p4->DeltaR(sublead_p4) < 0.3) sublead_is_prompt = true;
              }

            }

            //cout << lead_is_prompt << " " << sublead_is_prompt << endl;
            if (lead_is_prompt && sublead_is_prompt) {
              bkgCat = promptprompt;
            } else if (lead_is_prompt || sublead_is_prompt) {
              bkgCat = promptfake;
            } else {
	      bkgCat = fakefake;
	    }
          }
	  */

	  //Fill histograms according to diphoton or single photon category, as appropriate

          int dipho_category = l.DiphotonCategory(diphoton_index.first, diphoton_index.second, Higgs.Pt(), 2, 2, 2);
          int leadpho_category = l.PhotonCategory(diphoton_index.first, 2, 2);
          int subleadpho_category = l.PhotonCategory(diphoton_index.second, 2, 2);
/*
          if (bkgCat == promptfake) {
            l.FillHist("mass_pf",0, Higgs.M(), weight);
            l.FillHist("mass_pf",dipho_category+1, Higgs.M(), weight);
          } else if (bkgCat == fakefake) {
            l.FillHist("mass_ff",0, Higgs.M(), weight);
            l.FillHist("mass_ff",dipho_category+1, Higgs.M(), weight);
          } else {
            l.FillHist("mass",0, Higgs.M(), weight);
            l.FillHist("mass",dipho_category+1, Higgs.M(), weight);
          }
*/
	  //Only fill histograms for QCD and GJet if bkgCat is not promptprompt
	  //if ((l.itype[l.current]==1 || l.itype[l.current]==2 || l.itype[l.current]==3) && bkgCat==promptprompt) return;

          l.FillHist("pt",0, Higgs.Pt(), weight);
          l.FillHist("ptOverM",0, Higgs.Pt()/Higgs.M(), weight);
	  l.FillHist("eta",0, Higgs.Eta(), weight);
	  l.FillHist("decayAngle",0, decayAngle, weight);
	  l.FillHist("helicityAngle",0, helicityAngle, weight);
          l.FillHist("pho1_pt",0,lead_p4.Pt(), weight);
          l.FillHist("pho2_pt",0,sublead_p4.Pt(), weight);
          l.FillHist("pho1_ptOverM",0,lead_p4.Pt()/Higgs.M(), weight);
          l.FillHist("pho2_ptOverM",0,sublead_p4.Pt()/Higgs.M(), weight);
          l.FillHist("pho1_eta",0,lead_p4.Eta(), weight);
          l.FillHist("pho2_eta",0,sublead_p4.Eta(), weight);
          l.FillHist("pho_r9",0,l.pho_r9[diphoton_index.first], weight);
          l.FillHist("pho_r9",0,l.pho_r9[diphoton_index.second], weight);
          l.FillHist("maxeta",0,maxeta, weight);
          l.FillHist("ht",0, HT, weight);
          l.FillHist2D("ht_vs_m",0, HT, mH, weight);

          l.FillHist("pt",dipho_category+1, Higgs.Pt(), weight);
          l.FillHist("ptOverM",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
	  l.FillHist("eta",dipho_category+1, Higgs.Eta(), weight);
	  l.FillHist("decayAngle",dipho_category+1, decayAngle, weight);
	  l.FillHist("helicityAngle",dipho_category+1, helicityAngle, weight);
          l.FillHist("pho1_pt",dipho_category+1,lead_p4.Pt(), weight);
          l.FillHist("pho2_pt",dipho_category+1,sublead_p4.Pt(), weight);
          l.FillHist("pho1_ptOverM",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
          l.FillHist("pho2_ptOverM",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
          l.FillHist("pho1_eta",dipho_category+1,lead_p4.Eta(), weight);
          l.FillHist("pho2_eta",dipho_category+1,sublead_p4.Eta(), weight);
          l.FillHist("pho_r9",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
          l.FillHist("pho_r9",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
          l.FillHist("maxeta",dipho_category+1,maxeta, weight);
          l.FillHist("ht",dipho_category+1, HT, weight);
          l.FillHist2D("ht_vs_m",dipho_category+1, mH, HT, weight);

	  //Fill separately for low and high mass sidebands and for signal region

	  if (mH>100 && mH<110) {
	    l.FillHist("pt_mlow",0, Higgs.Pt(), weight);
	    l.FillHist("ptOverM_mlow",0, Higgs.Pt()/Higgs.M(), weight);
	    l.FillHist("eta_mlow",0, Higgs.Eta(), weight);
	    l.FillHist("decayAngle_mlow",0, decayAngle, weight);
	    l.FillHist("helicityAngle_mlow",0, helicityAngle, weight);
	    l.FillHist("pho1_pt_mlow",0,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_mlow",0,sublead_p4.Pt(), weight);
	    //l.FillHist("pho1_pt_mlow",0,lead_p4.Pt()*R_low, weight);
	    //l.FillHist("pho2_pt_mlow",0,sublead_p4.Pt()*R_low, weight);
	    l.FillHist("pho1_ptOverM_mlow",0,lead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho2_ptOverM_mlow",0,sublead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho1_eta_mlow",0,lead_p4.Eta(), weight);
	    l.FillHist("pho2_eta_mlow",0,sublead_p4.Eta(), weight);
	    l.FillHist("pho_r9_mlow",0,l.pho_r9[diphoton_index.first], weight);
	    l.FillHist("pho_r9_mlow",0,l.pho_r9[diphoton_index.second], weight);
	    l.FillHist("maxeta_mlow",0,maxeta, weight);
	    l.FillHist("ht_mlow",0, HT, weight);

	    l.FillHist("pt_mlow",dipho_category+1, Higgs.Pt(), weight);
	    l.FillHist("ptOverM_mlow",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
	    l.FillHist("eta_mlow",dipho_category+1, Higgs.Eta(), weight);
	    l.FillHist("decayAngle_mlow",dipho_category+1, decayAngle, weight);
	    l.FillHist("helicityAngle_mlow",dipho_category+1, helicityAngle, weight);
	    l.FillHist("pho1_pt_mlow",dipho_category+1,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_mlow",dipho_category+1,sublead_p4.Pt(), weight);
	    l.FillHist("pho1_ptOverM_mlow",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho2_ptOverM_mlow",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho1_eta_mlow",dipho_category+1,lead_p4.Eta(), weight);
	    l.FillHist("pho2_eta_mlow",dipho_category+1,sublead_p4.Eta(), weight);
	    l.FillHist("pho_r9_mlow",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
	    l.FillHist("pho_r9_mlow",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
	    l.FillHist("maxeta_mlow",dipho_category+1,maxeta, weight);
	    l.FillHist("ht_mlow",dipho_category+1, HT, weight);
	  } else if (mH>126 && mH<145) {
	    l.FillHist("pt_mhigh",0, Higgs.Pt(), weight);
	    l.FillHist("ptOverM_mhigh",0, Higgs.Pt()/Higgs.M(), weight);
	    l.FillHist("eta_mhigh",0, Higgs.Eta(), weight);
	    l.FillHist("decayAngle_mhigh",0, decayAngle, weight);
	    l.FillHist("helicityAngle_mhigh",0, helicityAngle, weight);
	    l.FillHist("pho1_pt_mhigh",0,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_mhigh",0,sublead_p4.Pt(), weight);
	    //l.FillHist("pho1_pt_mhigh",0,lead_p4.Pt()*R_high, weight);
	    //l.FillHist("pho2_pt_mhigh",0,sublead_p4.Pt()*R_high, weight);
	    l.FillHist("pho1_ptOverM_mhigh",0,lead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho2_ptOverM_mhigh",0,sublead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho1_eta_mhigh",0,lead_p4.Eta(), weight);
	    l.FillHist("pho2_eta_mhigh",0,sublead_p4.Eta(), weight);
	    l.FillHist("pho_r9_mhigh",0,l.pho_r9[diphoton_index.first], weight);
	    l.FillHist("pho_r9_mhigh",0,l.pho_r9[diphoton_index.second], weight);
	    l.FillHist("maxeta_mhigh",0,maxeta, weight);
	    l.FillHist("ht_mhigh",0, HT, weight);

	    l.FillHist("pt_mhigh",dipho_category+1, Higgs.Pt(), weight);
	    l.FillHist("ptOverM_mhigh",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
	    l.FillHist("eta_mhigh",dipho_category+1, Higgs.Eta(), weight);
	    l.FillHist("decayAngle_mhigh",dipho_category+1, decayAngle, weight);
	    l.FillHist("helicityAngle_mhigh",dipho_category+1, helicityAngle, weight);
	    l.FillHist("pho1_pt_mhigh",dipho_category+1,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_mhigh",dipho_category+1,sublead_p4.Pt(), weight);
	    l.FillHist("pho1_ptOverM_mhigh",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho2_ptOverM_mhigh",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho1_eta_mhigh",dipho_category+1,lead_p4.Eta(), weight);
	    l.FillHist("pho2_eta_mhigh",dipho_category+1,sublead_p4.Eta(), weight);
	    l.FillHist("pho_r9_mhigh",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
	    l.FillHist("pho_r9_mhigh",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
	    l.FillHist("maxeta_mhigh",dipho_category+1,maxeta, weight);
	    l.FillHist("ht_mhigh",dipho_category+1, HT, weight);
	  } else if (mH>=110 && mH<=126) {
	    l.FillHist("pt_msig",0, Higgs.Pt(), weight);
	    l.FillHist("ptOverM_msig",0, Higgs.Pt()/Higgs.M(), weight);
	    l.FillHist("eta_msig",0, Higgs.Eta(), weight);
	    l.FillHist("decayAngle_msig",0, decayAngle, weight);
	    l.FillHist("helicityAngle_msig",0, helicityAngle, weight);
	    l.FillHist("pho1_pt_msig",0,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_msig",0,sublead_p4.Pt(), weight);
	    l.FillHist("pho1_ptOverM_msig",0,lead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho2_ptOverM_msig",0,sublead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho1_eta_msig",0,lead_p4.Eta(), weight);
	    l.FillHist("pho2_eta_msig",0,sublead_p4.Eta(), weight);
	    l.FillHist("pho_r9_msig",0,l.pho_r9[diphoton_index.first], weight);
	    l.FillHist("pho_r9_msig",0,l.pho_r9[diphoton_index.second], weight);
	    l.FillHist("maxeta_msig",0,maxeta, weight);
	    l.FillHist("ht_msig",0, HT, weight);

	    l.FillHist("pt_msig",dipho_category+1, Higgs.Pt(), weight);
	    l.FillHist("ptOverM_msig",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
	    l.FillHist("eta_msig",dipho_category+1, Higgs.Eta(), weight);
	    l.FillHist("decayAngle_msig",dipho_category+1, decayAngle, weight);
	    l.FillHist("helicityAngle_msig",dipho_category+1, helicityAngle, weight);
	    l.FillHist("pho1_pt_msig",dipho_category+1,lead_p4.Pt(), weight);
	    l.FillHist("pho2_pt_msig",dipho_category+1,sublead_p4.Pt(), weight);
	    l.FillHist("pho1_ptOverM_msig",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho2_ptOverM_msig",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
	    l.FillHist("pho1_eta_msig",dipho_category+1,lead_p4.Eta(), weight);
	    l.FillHist("pho2_eta_msig",dipho_category+1,sublead_p4.Eta(), weight);
	    l.FillHist("pho_r9_msig",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
	    l.FillHist("pho_r9_msig",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
	    l.FillHist("maxeta_msig",dipho_category+1,maxeta, weight);
	    l.FillHist("ht_msig",dipho_category+1, HT, weight);
	  }
        }

	if(PADEBUG) 
		cout<<"myFillHistRed END"<<endl;
	
	
	if( runStatAnalysis ) {
		StatAnalysis(l,jentry);
	}

}


// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::StatAnalysis(LoopAll& l, Int_t jentry) 
{
}




// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
	vtxAna_.setBranchAdresses(t,"vtx_std_");
	vtxAna_.getBranches(t,"vtx_std_",s);
}


// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::PreselectPhotons(LoopAll& l, int jentry) 
{
	// Photon preselection
	pho_acc.clear();
	pho_presel.clear();
	pho_presel_lead.clear();
	pho_et.clear();
	l.pho_matchingConv->clear();

	// Nominal smearing
	corrected_pho_energy.clear(); corrected_pho_energy.resize(l.pho_n,0.); 
	int cur_type = l.itype[l.current];

	for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		std::vector<std::vector<bool> > p;
		PhotonReducedInfo phoInfo (
				*((TVector3*)l.pho_calopos->At(ipho)),
				((TLorentzVector*)l.pho_p4->At(ipho))->Energy(),
				energyCorrected[ipho],
				l.pho_isEB[ipho],
				l.pho_r9[ipho],
				false,
				(energyCorrectedError!=0?energyCorrectedError[ipho]:0)
					    );
		float pweight = 1.;
		float sweight = 1.;
		if( doEcorrectionSmear )  { 
		  eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.); 
		}
		if( cur_type == 0 ) {          // correct energy scale in data
		  eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
		  pweight *= sweight;
		}
		// apply mc-derived photon corrections, to data and MC alike
		corrected_pho_energy[ipho] = phoInfo.energy();
	}

	for(int ipho=0; ipho<l.pho_n; ++ipho) {

	  // match all photons in the original tree with the conversions from the merged collection and save the indices
	  int iConv  =l.matchPhotonToConversion(ipho);
	  if ( iConv>=0 )
		  (*l.pho_matchingConv).push_back(l.matchPhotonToConversion(ipho));
	  else
		  (*l.pho_matchingConv).push_back(-1);

	  // TLorentzVector * p4 = (TLorentzVector *) l.pho_p4->At(ipho);
	  TLorentzVector p4 = l.get_pho_p4(ipho,0,&corrected_pho_energy[0]);
	  // float eta  = fabs(((TVector3 *) l.pho_calopos->At(ipho))->Eta());
	  float eta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[ipho]))->Eta());
 	  // photon et wrt 0,0,0
	  float et = p4.Pt();
	  pho_et.push_back(et);
	  /// std::cerr << " " << p4->Pt() << " " << et << " " << eta;
	  
	  if( p4.Pt() < presel_scet2 || (eta>1.4442 && eta<1.566) || eta>presel_maxeta ) { 
		  /// std::cerr << std::endl;
	    continue;  
	  }
	  /// std::cerr << "keeping " << ipho << std::endl;
	  pho_acc.push_back(ipho);
	  
	  bool isEB = l.pho_isEB[ipho];
	  float & ecaliso = isEB ? presel_ecaliso_eb : presel_ecaliso_ee;
	  float & sieie = isEB ? presel_sieie_eb : presel_sieie_ee;
	  if( l.pho_ecalsumetconedr03[ipho] >= ecaliso ||  l.pho_sieie[ipho] >= sieie || l.pho_hoe[ipho] >= presel_hoe ) {
	    continue;
	  }
          
	  //FIXME trigger matching
	  pho_presel.push_back(ipho);
	} 

	std::sort(pho_acc.begin(),pho_acc.end(),
		  SimpleSorter<float,std::greater<float> >(&pho_et[0]));
	std::sort(pho_presel.begin(),pho_presel.end(),
		  SimpleSorter<float,std::greater<float> >(&pho_et[0]));

	if( pho_presel.size() > 1 ) {
		for(size_t ipho=0; ipho<pho_presel.size()-1; ++ipho ) {
			/// assert( ((TLorentzVector *)l.pho_p4->At(pho_presel[ipho]))->Pt() >= ((TLorentzVector *)l.pho_p4->At(pho_presel[ipho+1]))->Pt() );
			assert( pho_et[pho_presel[ipho]] >= pho_et[pho_presel[ipho+1]] );
		}
	}
	if( pho_acc.size()>1 ) {
		for(size_t ipho=0; ipho<pho_acc.size()-1; ++ipho ) {
			/// assert( ((TLorentzVector *)l.pho_p4->At(pho_acc[ipho]))->Pt() >= ((TLorentzVector *)l.pho_p4->At(pho_acc[ipho+1]))->Pt() );
			assert( pho_et[pho_acc[ipho]] >= pho_et[pho_acc[ipho+1]] );
		}
	}
	
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::FillReductionVariables(LoopAll& l, int jentry) 
{
	if(PADEBUG) 
		cout<<"myFillReduceVar START"<<endl;
  	
	PreselectPhotons(l,jentry);
		
	if(PADEBUG) 
		cout<<"myFillReduceVar END"<<endl;

}

// ----------------------------------------------------------------------------------------------------
bool PhotonAnalysis::SelectEventsReduction(LoopAll& l, int jentry) 
{

	if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
	// require at least two reconstructed photons to store the event
	if( pho_acc.size() < 2 || l.get_pho_p4( pho_acc[0], 0, &corrected_pho_energy[0] ).Pt() < presel_scet1 ) { return false; }
	
	///// int ipho1 = pho_acc[0];
	///// int ipho2 = pho_acc[1];
	///// 
	///// if( pho_presel.size() > 1 ) { 
	///// 	// use the first two preselected photons for the vertex algorithm
	///// 	ipho1 = pho_presel[0]; 
	///// 	ipho2 = pho_presel[1]; 
	///// } else if(pho_presel.size() > 0 ) {
	///// 	// if only one photon was preselected use the highest preselected and the higest non preselect photons 
	///// 	//    the event will be discarded by the analysis anyway
	///// 	ipho1 = pho_presel[0]; 
	///// 	ipho2 = pho_acc[0] == ipho1 ? pho_acc[1] : pho_acc[0]; 
	///// }
	///// assert( ipho1 != ipho2 );
	///// vtxAna_.clear();
	///// 
	///// if(PADEBUG)        cout << " SelectEventsReduction going to fill photon info " << endl;
	///// PhotonInfo pho1=l.fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions);
        ///// PhotonInfo pho2=l.fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions);
        ///// if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;
	///// 
        ///// // run vertex analysis
	///// l.vertexAnalysis(vtxAna_, pho1, pho2 );
        ///// // select vertex
	///// if( useDefaultVertex ) {
	///// 	l.vtx_std_ranked_list->clear();
	///// 	for(int ii=0;ii<l.vtx_std_n; ++ii) { l.vtx_std_ranked_list->push_back(ii); }
	///// } else {
	///// 	*l.vtx_std_ranked_list = l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames);
	///// 	if( l.vtx_std_ranked_list->size() != 0 ) {  
	///// 		l.vtx_std_sel = (*l.vtx_std_ranked_list)[0];
	///// 	} else {
	///// 		l.vtx_std_sel = 0;
	///// 		std::cerr << "NO VERTEX SELECTED " << l.event << " " << l.run << " " << std::endl;
	///// 	}
	///// 	// update the photons' pt
	///// 	for(int ipho=0; ipho<l.pho_n; ++ipho) {
	///// 		l.set_pho_p4(ipho, l.vtx_std_sel);
	///// 	}
	///// }
	
	vtxAna_.clear();
	l.vtx_std_ranked_list->clear();
	l.dipho_vtx_std_sel->clear();
	l.vtx_std_ranked_list->clear();
	l.vtx_std_sel=0;
	float maxSumPt = 0.;
	l.dipho_n = 0;

	// fill ID variables
	if( forcedRho >= 0. ) {
		l.rho = forcedRho;
	}
	l.FillCICInputs();
	l.FillCIC();

	if( pho_presel.size() < 2 ) {
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
		for(size_t id=0; id<diphotons.size(); ++id ) {
			
			int ipho1 = diphotons[id].first;
			int ipho2 = diphotons[id].second;
			
			if(PADEBUG)        cout << " SelectEventsReduction going to fill photon info " << endl;
			PhotonInfo pho1=l.fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions);
			PhotonInfo pho2=l.fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions);
			if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;
			
			l.vertexAnalysis(vtxAna_, pho1, pho2 );
			// make sure that vertex analysis indexes are in synch 
			assert( (int)id == vtxAna_.pairID(ipho1,ipho2) );
			
			l.vtx_std_ranked_list->push_back( l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames, mvaVertexSelection, 
									    tmvaPerVtxReader_, tmvaPerVtxMethod) );
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
			l.dipho_n = id+1;
			l.dipho_leadind[id] = diphotons[id].first;
			l.dipho_subleadind[id] = diphotons[id].second;
			l.dipho_vtxind[id] = l.dipho_vtx_std_sel->back();
			
			TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[id], l.dipho_vtxind[id], &corrected_pho_energy[0] );
			TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[id], l.dipho_vtxind[id], &corrected_pho_energy[0] );
			l.dipho_sumpt[id] = lead_p4.Pt() + sublead_p4.Pt();
			
			if( l.dipho_sumpt[id] > maxSumPt ) {
				l.vtx_std_sel = l.dipho_vtx_std_sel->back();
				maxSumPt = l.dipho_sumpt[id];
			}
		}
	}
	

	return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SkimEvents(LoopAll& l, int jentry)
{
	l.b_pho_n->GetEntry(jentry);
	if( l.pho_n < 2 ) {
		return false;
	}

	// do not run trigger selection on MC
	int filetype = l.itype[l.current];
	bool skipTrigger = !doTriggerSelection || filetype != 0 || triggerSelections.empty();
	if( ! skipTrigger ) {
		// get the trigger selection for this run 
		l.b_run->GetEntry(jentry);
		std::vector<TriggerSelection>::iterator isel = find(triggerSelections.begin(), triggerSelections.end(), l.run );
		if(isel == triggerSelections.end() ) {
			std::cerr << "No trigger selection for run " << l.run << "defined" << std::endl;
			return true;
		}

		// get the trigegr data
		l.b_hlt1_bit->GetEntry(jentry);
		l.b_hlt_path_names_HLT1->GetEntry(jentry);
		if( !  isel->pass(*(l.hlt_path_names_HLT1),*(l.hlt1_bit)) ) {
			/// std::cerr << "failed "  << std::endl;
			return false;
		}
		//// l.countersred[trigCounter_]++;
	}
	
	if( l.typerun == l.kReduce || l.typerun == l.kFillReduce ) {
		//// if( filetype == 2 ) { // photon+jet
		//// 	l.b_process_id->GetEntry(jentry);
		//// 	if( l.process_id == 18 ) {
		//// 		return false;
		//// 	}
		//// }
		
		if( filetype != 0 && ! (keepPP && keepPF && keepFF) ) {
			l.b_gp_n->GetEntry(jentry);
			l.b_gp_mother->GetEntry(jentry);
			l.b_gp_status->GetEntry(jentry);
			l.b_gp_pdgid->GetEntry(jentry);
			l.b_gp_p4->GetEntry(jentry);
			
			int np = 0;
			for(int ip=0;ip<l.gp_n;++ip) {
				if( l.gp_status[ip] != 1 || l.gp_pdgid[ip] != 22 ) { 
					continue; 
				}
				TLorentzVector * p4 = (TLorentzVector*) l.gp_p4->At(ip);
				if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
				int mother_id = abs( l.gp_pdgid[ l.gp_mother[ip] ] );
				if( mother_id <= 25 ) { ++np; }
				if( np >= 2 ) { break; }
			}
			if( np >= 2 && ! keepPP ) { return false; }
			if( np == 1 && ! keepPF ) { return false; }
			if( np == 0 && ! keepFF ) { return false; }
		}
	}

	return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
	return true;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
	vtxAna_.branches(outputTree,"vtx_std_");	

	l.pho_matchingConv = new  std::vector<int>();
	l.Branch_pho_matchingConv(outputTree);
	
	l.vtx_std_evt_mva = new std::vector<float>();
	l.vtx_std_ranked_list = new std::vector<std::vector<int> >();
	l.pho_tkiso_recvtx_030_002_0000_10_01 = new std::vector<std::vector<float> >();
	l.pho_cic6cutlevel_lead = new std::vector<std::vector<Short_t> >();
	l.pho_cic6passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
	l.pho_cic6cutlevel_sublead = new std::vector<std::vector<Short_t> >();
	l.pho_cic6passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
	l.pho_cic4cutlevel_lead = new std::vector<std::vector<Short_t> >();
	l.pho_cic4passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
	l.pho_cic4cutlevel_sublead = new std::vector<std::vector<Short_t> >();
	l.pho_cic4passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
	l.dipho_vtx_std_sel =  new std::vector<int>();

	l.Branch_vtx_std_evt_mva(outputTree);
	l.Branch_vtx_std_ranked_list(outputTree);
	l.Branch_vtx_std_sel(outputTree);
	l.Branch_pho_tkiso_recvtx_030_002_0000_10_01(outputTree);
	l.Branch_pho_tkiso_badvtx_040_002_0000_10_01(outputTree);
	l.Branch_pho_tkiso_badvtx_id(outputTree);
	l.Branch_pho_drtotk_25_99(outputTree);

	l.Branch_dipho_n(outputTree);
	l.Branch_dipho_leadind(outputTree);
	l.Branch_dipho_subleadind(outputTree);
	l.Branch_dipho_vtxind(outputTree);
	l.Branch_dipho_sumpt(outputTree);
	/// l.Branch_dipho_leadet(outputTree);
	/// l.Branch_dipho_subleadet(outputTree);
	/// l.Branch_dipho_leadeta(outputTree);
	/// l.Branch_dipho_subleadeta(outputTree);
	/// l.Branch_dipho_leadci6cindex(outputTree);
	/// l.Branch_dipho_subleadci6cindex(outputTree);
	/// l.Branch_dipho_leadci4cindex(outputTree);
	/// l.Branch_dipho_subleadci4cindex(outputTree);
	/// l.Branch_dipho_mass(outputTree);
	/// l.Branch_dipho_pt(outputTree);
	/// l.Branch_dipho_eta(outputTree);
	/// l.Branch_dipho_phi(outputTree);
	/// l.Branch_dipho_cts(outputTree);
	
	l.Branch_pho_cic6cutlevel_lead( outputTree );
	l.Branch_pho_cic6passcuts_lead( outputTree );
	l.Branch_pho_cic6cutlevel_sublead( outputTree );
	l.Branch_pho_cic6passcuts_sublead( outputTree );
	l.Branch_pho_cic4cutlevel_lead( outputTree );
	l.Branch_pho_cic4passcuts_lead( outputTree );
	l.Branch_pho_cic4cutlevel_sublead( outputTree );
	l.Branch_pho_cic4passcuts_sublead( outputTree );
	
}

