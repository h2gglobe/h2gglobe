#include "../interface/CategoryAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
CategoryAnalysis::CategoryAnalysis()  : 
    name_("CategoryAnalysis")
{
}

// ----------------------------------------------------------------------------------------------------
CategoryAnalysis::~CategoryAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void CategoryAnalysis::Term(LoopAll& l) 
{

    if (! l.is_subjob){ // no need to waste time when running a subjob
        std::string outputfilename = (std::string) l.histFileName;
    }

    eventListText.close();
    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;

}

// ----------------------------------------------------------------------------------------------------
void CategoryAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
        cout << "InitRealCategoryAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;

    //  std::string outputfilename = (std::string) l.histFileName;
    eventListText.open(Form("%s",l.outputTextFileName.c_str()));
    FillSignalLabelMap(l);
    //
    // These parameters are set in the configuration file
    std::cout
        << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << "CategoryAnalysis " << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << "leadEtCut "<< leadEtCut << "\n"
        << "subleadEtCut "<< subleadEtCut << "\n"
        << "doTriggerSelection "<< doTriggerSelection << "\n"
        << "nEtaCategories "<< nEtaCategories << "\n"
        << "nR9Categories "<< nR9Categories << "\n"    
        << "nPtCategories "<< nPtCategories << "\n"  
        << "efficiencyFile " << efficiencyFile << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << std::endl;

    PhotonAnalysis::Init(l);
    PhotonAnalysis::doSystematics = false;

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
        l.histoContainer[ind].setScale(1.);
    }

    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    nInclusiveCategories_ = 4;

  
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
        eResolSmearer->doEnergy(false); // allows for future reweighting also
        eResolSmearer->scaleOrSmear(false);
        photonSmearers_.push_back(eResolSmearer);
    }
    if( doRegressionSmear ) {
        // energy regression. smearing
        std::cerr << __LINE__ << std::endl; 
        eRegressionSmearer = new EnergySmearer( eSmearPars );
        eRegressionSmearer->name("regSig");
        eRegressionSmearer->doEnergy(false);// allows for future reweighting also
        eRegressionSmearer->doRegressionSigma(true);
        photonSmearers_.push_back(eRegressionSmearer);
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
    if( doPhotonMvaIdSmear ) {
        // trigger efficiency
        std::cerr << __LINE__ << std::endl; 
        photonMvaIdSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
        photonMvaIdSmearer->name("phoIdMva");
        photonMvaIdSmearer->setEffName("effL1HLT");
        photonMvaIdSmearer->doVtxEff(false);
        photonMvaIdSmearer->doMvaIdEff(true);
        photonMvaIdSmearer->init();
        diPhotonSmearers_.push_back(photonMvaIdSmearer);
    }
    if(doKFactorSmear) {
        // kFactor efficiency
        std::cerr << __LINE__ << std::endl; 
        kFactorSmearer = new KFactorSmearer( kfacHist );
        kFactorSmearer->name("kFactor");
        kFactorSmearer->init();
        genLevelSmearers_.push_back(kFactorSmearer);
    }
    if(doInterferenceSmear) {
        // interference efficiency
        std::cerr << __LINE__ << std::endl; 
        interferenceSmearer = new InterferenceSmearer(2.5e-2,0.);
        genLevelSmearers_.push_back(interferenceSmearer);
    }

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed
    int nVBFCategories   = ((int)includeVBF)*nVBFEtaCategories*nVBFDijetJetCategories;
    std::sort(bdtCategoryBoundaries.begin(),bdtCategoryBoundaries.end(), std::greater<float>() );

    // Make sure the Map is filled
    FillSignalLabelMap(l);

    // Initialize all MVA ---------------------------------------------------//
    l.SetAllMVA();
    // UCSD
    l.tmvaReaderID_UCSD->BookMVA("Gradient"      ,photonLevelMvaUCSD.c_str()  );
    l.tmvaReader_dipho_UCSD->BookMVA("Gradient"  ,eventLevelMvaUCSD.c_str()   );
    // New ID MVA
    l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevelNewIDMVA_EB.c_str());
    l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevelNewIDMVA_EE.c_str());
    // MIT 
    l.tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost",photonLevelMvaMIT_EB.c_str());
    l.tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost",photonLevelMvaMIT_EE.c_str());
    l.tmvaReader_dipho_MIT->BookMVA("Gradient"   ,eventLevelMvaMIT.c_str()    );
    // ----------------------------------------------------------------------//

    if(PADEBUG) 
        cout << "InitRealCategoryAnalysis END"<<endl;

    // FIXME book of additional variables
}

bool CategoryAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4,
				float & mass, float & evweight, int & category, int & diphoton_id, bool & isCorrectVertex, float &kinematic_bdtout,
				bool isSyst, 
				float syst_shift, bool skipSelection, 
				BaseGenLevelSmearer *genSys, BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys) 
{

    assert( isSyst || ! skipSelection );

    int cur_type = l.itype[l.current];
    float sampleweight = l.sampleContainer[l.current_sample_index].weight;
    /// diphoton_id = -1;
    
    std::pair<int,int> diphoton_index;
    std::pair<int,int> diphoton_index_cic;
    std::pair<int,int> diphoton_index_kinonly;
   
    // do gen-level dependent first (e.g. k-factor); only for signal
    genLevWeight=1.;
    if(cur_type!=0 ) {
	applyGenLevelSmearings(genLevWeight,gP4,l.pu_n,cur_type,genSys,syst_shift);
    }
	
    // first apply corrections and smearing on the single photons 
    smeared_pho_energy.clear(); smeared_pho_energy.resize(l.pho_n,0.); 
    smeared_pho_r9.clear();     smeared_pho_r9.resize(l.pho_n,0.); 
    smeared_pho_weight.clear(); smeared_pho_weight.resize(l.pho_n,1.);
    applySinglePhotonSmearings(smeared_pho_energy, smeared_pho_r9, smeared_pho_weight, cur_type, l, energyCorrected, energyCorrectedError,
			       phoSys, syst_shift);

    // Exclusive Modes
    int diphotonVBF_id = -1;
    int ijet1, ijet2;
    VBFevent = false;
    int diphotonVBF_id_CIC=-1;
    bool VBFevent_CIC = false;
    double subleadEtVBFCut_CIC=25.;

    // VBF
    if((includeVBF || includeVHhad)&&l.jet_algoPF1_n>1 && !isSyst /*avoid rescale > once*/) {
      l.RescaleJetEnergy();
    }

    // Flag whether this is a VBF event (separately for MVA and CIC because sublead Et cut is different)
    if(includeVBF) {
      diphotonVBF_id = l.DiphotonMITPreSelection(leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] );
      VBFevent=VBFTag2012(ijet1, ijet2, l, diphotonVBF_id, &smeared_pho_energy[0], false, 0., 0.);

      diphotonVBF_id_CIC = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut_CIC, 4,false, &smeared_pho_energy[0], true); 
      VBFevent_CIC=VBFTag2012(ijet1, ijet2, l, diphotonVBF_id_CIC, &smeared_pho_energy[0], false, 0., 0.);
    }


    // Determine whether or not event passes CIC selection in each category

    bool passCiC = false;
    int catCiC = -1;

    int diphoton_id_cic = -1;
    diphoton_id_cic = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, nPhotonCategories_,applyPtoverM, &smeared_pho_energy[0] ); 

    if (diphoton_id_cic > -1 && !VBFevent_CIC) {

      diphoton_index_cic = std::make_pair( l.dipho_leadind[diphoton_id_cic],  l.dipho_subleadind[diphoton_id_cic] );
      // bring all the weights together: lumi & Xsection, single gammas, pt kfactor
      float evweight_cic = weight * smeared_pho_weight[diphoton_index_cic.first] * smeared_pho_weight[diphoton_index_cic.second] * genLevWeight;

      TLorentzVector lead_p4_cic = l.get_pho_p4( l.dipho_leadind[diphoton_id_cic], l.dipho_vtxind[diphoton_id_cic], &smeared_pho_energy[0]);
      TLorentzVector sublead_p4_cic = l.get_pho_p4( l.dipho_subleadind[diphoton_id_cic], l.dipho_vtxind[diphoton_id_cic], &smeared_pho_energy[0]);
      TLorentzVector Higgs_cic = lead_p4_cic + sublead_p4_cic;   

      int selectioncategory_cic = l.DiphotonCategory(diphoton_index_cic.first,diphoton_index_cic.second,Higgs_cic.Pt(),nEtaCategories,nR9Categories,0);
      float mass_cic    = Higgs_cic.M();

      assert( evweight_cic >= 0. ); 

      l.FillHist("mass_cic",0, mass_cic, evweight_cic);

      if( mass_cic>=massMin && mass_cic<=massMax ) {
	passCiC = true;
	catCiC = selectioncategory_cic;
      }

    }


    // Determine whether or not event passes MVA cut in each category with two additional categories:
    // pass kinematic preselection (pt and eta cuts only) but fail preselection,
    // and pass preselection but fail MVA cut

    bool passKin = false;
    int selectioncategory_kinonly = -1;
    bool passMVApresel = false;
    bool passMVA = false;
    int selectioncategory = -1;
    int catMVA = -2;
    float diphobdt_output = -2.;

    // Kinematic preselection (pt and eta cuts only)
    int diphoton_id_kinonly = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0], true );
    if (diphoton_id_kinonly > -1 ) {
      passKin = true;
      diphoton_index_kinonly = std::make_pair( l.dipho_leadind[diphoton_id_kinonly],  l.dipho_subleadind[diphoton_id_kinonly] );
      TLorentzVector lead_p4_kinonly = l.get_pho_p4( l.dipho_leadind[diphoton_id_kinonly], l.dipho_vtxind[diphoton_id_kinonly], &smeared_pho_energy[0]);
      TLorentzVector sublead_p4_kinonly = l.get_pho_p4( l.dipho_subleadind[diphoton_id_kinonly], l.dipho_vtxind[diphoton_id_kinonly], &smeared_pho_energy[0]);
      TLorentzVector Higgs_kinonly = lead_p4_kinonly + sublead_p4_kinonly;   
      selectioncategory_kinonly = l.DiphotonCategory(diphoton_index_kinonly.first,diphoton_index_kinonly.second,Higgs_kinonly.Pt(),nEtaCategories,nR9Categories,0);
    }

    // Preselection
    diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] ); 
    if (diphoton_id > -1 ) {
      passMVApresel = true;
      catMVA = -1;
    } else {
      diphoton_id = diphoton_id_cic;
    }

    // Calculate diphoon MVA output for events which pass either MVA preselection or CIC selection
    if (diphoton_id > -1 ) {

      diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );

      // bring all the weights together: lumi & Xsection, single gammas, pt kfactor
      evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
      if( ! isSyst ) {
	l.countersred[diPhoCounter_]++;
      }

      TLorentzVector lead_p4, sublead_p4, Higgs;
      float lead_r9, sublead_r9;
      TVector3 * vtx;
      fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  

      float mass    = Higgs.M();
      float ptHiggs = Higgs.Pt();

      // Mass Resolution of the Event
      // Mass Resolution of the Event
      ///// double beamspotSigma=-100;
      ///// if(l.version<13) {
      ///// 	  beamspotSigma=5.8;
      ///// } else {
      ///// 	  beamspotSigma=4.8;
      ///// }
      ///// massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories,beamspotSigma);
      massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories,beamspotWidth);
      float vtx_mva  = l.vtx_std_evt_mva->at(diphoton_id);
      float sigmaMrv = massResolutionCalculator->massResolutionEonly();
      float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
      float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
      // easy to calculate vertex probability from vtx mva output
      float vtxProb   = 1.-0.49*(vtx_mva+1.0); /// should better use this: vtxAna_.setPairID(diphoton_id); vtxAna_.vertexProbability(vtx_mva); PM

      float phoid_mvaout_lead = l.photonIDMVANew(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4,bdtTrainingPhilosophy.c_str()) + photonIDMVAShift_EB;
      float phoid_mvaout_sublead = l.photonIDMVANew(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4,bdtTrainingPhilosophy.c_str()) + photonIDMVAShift_EE;

      // apply di-photon level smearings and corrections
      selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
      if( cur_type != 0 && doMCSmearing ) {
	applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, phoid_mvaout_lead,phoid_mvaout_sublead,
			       diPhoSys, syst_shift);
	isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
      }
                           
      // Must be calculated after photon id has potentially been smeared
      diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id] ,vtxProb,lead_p4,sublead_p4 ,sigmaMrv,sigmaMwv,sigmaMeonly ,bdtTrainingPhilosophy.c_str() ,phoid_mvaout_lead,phoid_mvaout_sublead);
      kinematic_bdtout = diphobdt_output;

      bool isEBEB  = (lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);
      int category = -2;
      if (passMVApresel) category = GetBDTBoundaryCategory(diphobdt_output, isEBEB,VBFevent);

      assert( evweight >= 0. ); 

      if( ! isSyst ) {

	l.FillCounter( "Accepted", weight );
	l.FillCounter( "Smeared", evweight );
	sumaccept += weight;
	sumsmear += evweight;

	if (!VBFevent) {

	  if (passMVApresel) {
	    if (diphobdt_output>-0.05) {
	      l.FillHist("mass_mva",0, mass, evweight);
	    } else {
	      l.FillHist("mass_mva_bdtoutlt05",0, mass, evweight);
	    }
	  } else if (passCiC) {
	    if (diphobdt_output>-0.05) {
	      l.FillHist("mass_failpresel_bdtoutgt05",0, mass, evweight);
	    } else {
	      l.FillHist("mass_failpresel_bdtoutlt05",0, mass, evweight);
	    }
	  }

	  if( mass>=massMin && mass<=massMax) {

	    if (diphobdt_output > -0.05) {
	      float deltaEta = lead_p4.Eta()-sublead_p4.Eta();     
	      float max_eta = max(fabs(lead_p4.Eta()),fabs(sublead_p4.Eta()));
	      float min_r9  = min(lead_r9,sublead_r9);
	      l.FillHist2D("eta1_vs_deltaEta",category+2,fabs(deltaEta),fabs(lead_p4.Eta()),evweight);
	      l.FillHist2D("eta2_vs_eta1",category+2,lead_p4.Eta(),sublead_p4.Eta(),evweight);
	      l.FillHist2D("minR9_vs_maxEta",category+2,max_eta,min_r9,evweight);
	    }

	    if (passMVApresel) {
	      l.FillHist2D("pt_vs_bdtout",selectioncategory,diphobdt_output,ptHiggs,evweight);
	      l.FillHist("bdtout",selectioncategory,diphobdt_output,evweight);
	      if (ptHiggs<40.) {
		l.FillHist("bdtout_lowPt",selectioncategory,diphobdt_output,evweight);
	      } else {
		l.FillHist("bdtout_highPt",selectioncategory,diphobdt_output,evweight);
	      }

	      if (passCiC) {
		l.FillHist("bdtout_passCiC",selectioncategory,diphobdt_output,evweight);
		if (ptHiggs<40.) {
		  l.FillHist("bdtout_passCiC_lowPt",selectioncategory,diphobdt_output,evweight);
		} else {
		  l.FillHist("bdtout_passCiC_highPt",selectioncategory,diphobdt_output,evweight);
		}
	      }
	    }

	    if (category>-1 && category<4) {
	      passMVA = true;
	      catMVA = category;
	    } else if (catCiC==0) {
	      eventListText <<  " Run=" << l.run << "  LS=" << l.lumis << "  Event=" << l.event << " BDTCAT=" << category << " ggM=" << mass << " gg_Pt=" << ptHiggs << " LeadPhotonPhoid=" << phoid_mvaout_lead << " SubleadPhotonPhoid=" <<phoid_mvaout_sublead << " diphotonBDT=" << diphobdt_output << " photon1Eta=" << lead_p4.Eta() <<" photon2Eta="<<sublead_p4.Eta() << " sigmaMrv="<<sigmaMrv << " sigmaMwv=" << sigmaMwv << " photon1Pt="<<lead_p4.Pt()<<" photon2Pt="<<sublead_p4.Pt() << " vtxProb="<<vtxProb <<" cosDphi="<<TMath::Cos(lead_p4.Phi() - sublead_p4.Phi()) << " r9_1=" <<lead_r9 <<" r9_2=" <<sublead_r9  <<" E1="<<lead_p4.E()<<" E2="<<sublead_p4.E() << endl;
	    }

	  }

	}

      }

    } else if (catCiC==0) {
      eventListText <<  " Run=" << l.run << "  LS=" << l.lumis << "  Event=" << l.event << endl;
    }

    if (!VBFevent) {
      int diphocat = selectioncategory;
      if (selectioncategory == -1) diphocat = selectioncategory_kinonly;
      int passLevelMVA = 0;
      if (passKin) passLevelMVA = 1;
      if (passMVApresel) passLevelMVA = 2;
      if (passMVA) passLevelMVA = 3;
      l.FillHist2D("diphocat_vs_passMVA",0, float(passLevelMVA), float(diphocat), weight*genLevWeight);
    }

    l.FillHist2D("passCiC_vs_passMVA",0, float(passMVA), float(passCiC), weight*genLevWeight);
    l.FillHist2D("passCiC_vs_passMVA_cat",0, float(catMVA), float(catCiC), weight*genLevWeight);

    bool VBFevent_MVA = false;
    if (VBFevent && diphobdt_output>-0.05) VBFevent_MVA = true;
    l.FillHist2D("passVBFCiC_vs_passVBFMVA",0, float(VBFevent_MVA), float(VBFevent_CIC), weight*genLevWeight);

    return false;
}


int CategoryAnalysis::category(std::vector<float> & v, float val)
{
	if( val == v[0] ) { return 0; }
	std::vector<float>::iterator bound =  lower_bound( v.begin(), v.end(), val, std::greater<float>  ());
	int cat = ( val >= *bound ? bound - v.begin() - 1 : bound - v.begin() );
    //for (int i=0; i<v.size(); i++) cout << v[i] << endl;
	if( cat >= v.size() - 1 ) { cat = -1; }
	return cat;
}


// ----------------------------------------------------------------------------------------------------
int CategoryAnalysis::GetBDTBoundaryCategory(float bdtout, bool isEB, bool VBFevent){

    if (bdtTrainingPhilosophy=="UCSD"){
        if (isEB) { // 6 Categories for the EB-EB 
            if (bdtout < -0.30) return 5;
            if (bdtout >=-0.30 && bdtout < 0.00) return 4;
            if (bdtout >= 0.00 && bdtout < 0.30) return 3;
            if (bdtout >= 0.30 && bdtout < 0.60) return 2;
            if (bdtout >= 0.60 && bdtout < 0.70) return 1;
            if (bdtout >= 0.70) return 0;
        }
        else {// 2 Categories for the EB/EE 
            if (bdtout <  0.1) return 7;
            if (bdtout >= 0.1) return 6;
        }

    } else if (bdtTrainingPhilosophy=="MIT"){
	int cat = category( bdtCategoryBoundaries, bdtout );
	if( VBFevent && cat > -1 ) cat = bdtCategoryBoundaries.size();
	return cat;
    } else std::cerr << "No BDT Philosophy known - " << bdtTrainingPhilosophy << std::endl;
}


void CategoryAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    eResolSmearer->resetRandom();
}



// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
