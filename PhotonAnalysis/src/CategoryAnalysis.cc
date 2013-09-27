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
	<< "efficiencyFileMVA " << efficiencyFileMVA << "\n"
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

    effSmearPars.categoryType = effPhotonCategoryType;
    effSmearPars.n_categories = effPhotonNCat;
    effSmearPars.efficiency_file = efficiencyFileMVA;

    diPhoEffSmearPars.n_categories = 8;
    diPhoEffSmearPars.efficiency_file = efficiencyFileMVA;

    if( doEcorrectionSmear ) {
        // instance of this smearer done in PhotonAnalysis
        photonSmearers_.push_back(eCorrSmearer);
    }
    if( doEscaleSmear ) {
    setupEscaleSmearer();
    }
    if( doEresolSmear ) {
        // energy resolution smearing
    setupEresolSmearer();
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
    if(doPdfWeightSmear) {
        // PdfWeights efficiency (For now only consider QCD Scale Uncertainty 
        std::cerr << __LINE__ << std::endl; 
        pdfWeightSmearer = new PdfWeightSmearer( pdfWeightHist,"up","down");
        pdfWeightSmearer->name("pdfWeight");
        pdfWeightSmearer->init();
        genLevelSmearers_.push_back(pdfWeightSmearer);
    }
    if(doInterferenceSmear) {
        // interference efficiency
        std::cerr << __LINE__ << std::endl; 
        interferenceSmearer = new InterferenceSmearer(2.5e-2,0.);
        genLevelSmearers_.push_back(interferenceSmearer);
    }

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed
    // initialize the analysis variables
    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;

    std::sort(bdtCategoryBoundaries.begin(),bdtCategoryBoundaries.end(), std::greater<float>() );
    nInclusiveCategories_ = bdtCategoryBoundaries.size()-1;

    nVBFCategories   = ((int)includeVBF)*( (mvaVbfSelection && !multiclassVbfSelection) ? mvaVbfCatBoundaries.size()-1 : nVBFEtaCategories*nVBFDijetJetCategories );
    if(includeVHlep){
        nVHlepCategories = nElectronCategories + nMuonCategories;
    }
    if(includeVHmet){
    nVHmetCategories = nMetCategories;
    }
    
    std::sort(mvaVbfCatBoundaries.begin(),mvaVbfCatBoundaries.end(), std::greater<float>() );
    if (multiclassVbfSelection) {
        std::vector<int> vsize;
        vsize.push_back((int)multiclassVbfCatBoundaries0.size());
        vsize.push_back((int)multiclassVbfCatBoundaries1.size());
        vsize.push_back((int)multiclassVbfCatBoundaries2.size());
        std::sort(vsize.begin(),vsize.end(), std::greater<int>());
        // sanity check: there sould be at least 2 vectors with size==2
        if (vsize[0]<2 || vsize[1]<2 ){
            std::cout << "Not enough category boundaries:" << std::endl;
            std::cout << "multiclassVbfCatBoundaries0 size = " << multiclassVbfCatBoundaries0.size() << endl;
            std::cout << "multiclassVbfCatBoundaries1 size = " << multiclassVbfCatBoundaries1.size() << endl;
            std::cout << "multiclassVbfCatBoundaries2 size = " << multiclassVbfCatBoundaries2.size() << endl;
            assert( 0 );
        }
        nVBFCategories   = vsize[0]-1;
        cout << "@@@@@@@@@@@@@@@@@ 	nVBFCategories = " << 	nVBFCategories << endl;
        std::sort(multiclassVbfCatBoundaries0.begin(),multiclassVbfCatBoundaries0.end(), std::greater<float>() );
        std::sort(multiclassVbfCatBoundaries1.begin(),multiclassVbfCatBoundaries1.end(), std::greater<float>() );
        std::sort(multiclassVbfCatBoundaries2.begin(),multiclassVbfCatBoundaries2.end(), std::greater<float>() );
    }

    nCategories_=(nInclusiveCategories_+nVBFCategories+nVHlepCategories+nVHmetCategories);


    // Make sure the Map is filled
    FillSignalLabelMap(l);

    // Initialize all MVA ---------------------------------------------------//
    l.SetAllMVA();
    // UCSD
    l.tmvaReaderID_UCSD->BookMVA("Gradient"      ,photonLevelMvaUCSD.c_str()  );
    l.tmvaReader_dipho_UCSD->BookMVA("Gradient"  ,eventLevelMvaUCSD.c_str()   );
    // New ID MVA
    if( photonLevelNewIDMVA_EB != "" && photonLevelNewIDMVA_EE != "" ) {
    l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevelNewIDMVA_EB.c_str());
    l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevelNewIDMVA_EE.c_str());
    } else { 
    assert( dataIs2011 );
    }
    // MIT 
    if( photonLevelMvaMIT_EB != "" && photonLevelMvaMIT_EE != "" ) {
    l.tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost",photonLevelMvaMIT_EB.c_str());
    l.tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost",photonLevelMvaMIT_EE.c_str());
    } else {
    assert( ! dataIs2011 );
    }
    l.tmvaReader_dipho_MIT->BookMVA("Gradient"   ,eventLevelMvaMIT.c_str()    );
    // ----------------------------------------------------------------------//

    if(PADEBUG) 
        cout << "InitRealCategoryAnalysis END"<<endl;

    cout << "------- BUTTERY BUS ---- " << endl;
    cout << "reweighBS - " << reweighBeamspot << endl;
    cout << "saveDatTr - " << saveDatacardTrees_ << endl;
    cout << "-------------------------" << endl;
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
    float sampleweight = l.sampleContainer[l.current_sample_index].weight();
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

    int mu_ind=-1;
    int el_ind=-1;

    int muVtx=-1;
    int elVtx=-1;
    
    // Exclusive Modes
    int diphotonVBF_id = -1;
    vbfIjet1=-1, vbfIjet2=-1;
    VBFevent = false;
        
    int diphotonVBF_id_CIC=-1;
    bool VBFevent_CIC = false;

    int diphotonVHlep_id = -1;
    VHmuevent = false;
    VHmuevent_cat=0;
    VHelevent = false;
    VHelevent_cat=0;
        
    int diphotonVHmet_id = -1; 
    VHmetevent = false; 
    VHmetevent_cat=0; 

    int diphotonVHhad_id = -1;
    VHhadevent = false;

    if(includeVHlep){
	VHmuevent=MuonTag2012B(l, diphotonVHlep_id, mu_ind, muVtx, VHmuevent_cat, &smeared_pho_energy[0], lep_sync, true, phoidMvaCut, 0., smeared_pho_weight, false);
	ElectronStudies2012B(l, &smeared_pho_energy[0], true,  phoidMvaCut, 0., 0., jentry);
	int diphotonVH_ele_id=-1;
	VHelevent=ElectronTag2012B(l, diphotonVH_ele_id, el_ind, elVtx, VHelevent_cat, &smeared_pho_energy[0], lep_sync, true, phoidMvaCut, 0., smeared_pho_weight, false);
	if(!VHmuevent && VHelevent){
	    diphotonVHlep_id=diphotonVH_ele_id;
	}
    }

    if(includeVHmet && !dataIs2011) {
	//	    std::cout << "+++PFMET UNCORR " << l.met_pfmet << std::endl;
	if(!isSyst) VHmetevent=METTag2012B(l, diphotonVHmet_id, VHmetevent_cat, &smeared_pho_energy[0], met_sync, true, phoidMvaCut, false); 
	if(isSyst)  VHmetevent=METTag2012B(l, diphotonVHmet_id, VHmetevent_cat, &smeared_pho_energy[0], met_sync, true, phoidMvaCut, true); 
	// FIXME  need to un-flag events failing the diphoton mva cut.
    }

    // VBF
    if((includeVBF || includeVHhad)&&l.jet_algoPF1_n>1 && !isSyst /*avoid rescale > once*/) {
	l.RescaleJetEnergy();
    }

    // Flag whether this is a VBF event (separately for MVA and CIC because sublead Et cut is different)
    if(includeVBF) {   
	diphotonVBF_id = l.DiphotonMITPreSelection(leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0] );
	VBFevent= VBFTag2012(vbfIjet1, vbfIjet2, l, diphotonVBF_id, &smeared_pho_energy[0], false, 0., 0.);

	diphotonVBF_id_CIC = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVBFCut, subleadEtVBFCut, 4,false, &smeared_pho_energy[0], true); 
	VBFevent_CIC=VBFTag2012(vbfIjet1, vbfIjet2, l, diphotonVBF_id_CIC, &smeared_pho_energy[0], false, 0., 0.);
    }

    //Exclude events that are classified as lepton tag or MET tag
    bool isOtherExclusive=false;
    if(includeVHlep&&VHmuevent){
	isOtherExclusive=true;
    } else if (includeVHlep&&VHelevent){
	isOtherExclusive=true;
    } else if(includeVBF&&VBFevent) {
	isOtherExclusive=false;
    } else if(includeVHmet&&VHmetevent) {
	isOtherExclusive=true;
    }
    if (isOtherExclusive) return false;

    // Determine whether or not event passes CIC selection in each category

    bool passCiC = false;
    int catCiC = -1;

    int diphoton_id_cic = -1;
    std::vector<bool> veto_indices;
    veto_indices.clear();
    diphoton_id_cic = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0], false, -1, veto_indices, cicCutLevels );

    if (diphoton_id_cic > -1 && !VBFevent_CIC) {

      diphoton_index_cic = std::make_pair( l.dipho_leadind[diphoton_id_cic],  l.dipho_subleadind[diphoton_id_cic] );
      // bring all the weights together: lumi & Xsection, single gammas, pt kfactor
      float evweight_cic = weight * smeared_pho_weight[diphoton_index_cic.first] * smeared_pho_weight[diphoton_index_cic.second] * genLevWeight;

      TLorentzVector lead_p4_cic, sublead_p4_cic, Higgs_cic;
      float lead_r9_cic, sublead_r9_cic;
      TVector3 * vtx_cic;
      fillDiphoton(lead_p4_cic, sublead_p4_cic, Higgs_cic, lead_r9_cic, sublead_r9_cic, vtx_cic, &smeared_pho_energy[0], l, diphoton_id_cic);
      int selectioncategory_cic = l.DiphotonCategory(diphoton_index_cic.first,diphoton_index_cic.second,Higgs_cic.Pt(),nEtaCategories,nR9Categories,0);

      if( cur_type != 0 && doMCSmearing ) {
	  applyDiPhotonSmearings(Higgs_cic, *vtx_cic, selectioncategory_cic, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight_cic, zero_, zero_,diPhoSys, syst_shift);
      }

      float mass_cic    = Higgs_cic.M();

      // apply beamspot reweighting if necessary
      if(reweighBeamspot && cur_type!=0) {
	  evweight_cic*=BeamspotReweight(vtx_cic->Z(),((TVector3*)l.gv_pos->At(0))->Z());
      }

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
    mass = 0.;

    // Kinematic preselection (pt and eta cuts only)
    int diphoton_id_kinonly = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0], false, true, -100, -1, false);
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
        
      // apply beamspot reweighting if necessary
      if(reweighBeamspot && cur_type!=0) {
	  evweight*=BeamspotReweight(vtx->Z(),((TVector3*)l.gv_pos->At(0))->Z());
      }

      massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,massResoPars,nR9Categories,nEtaCategories,beamspotSigma);
      float vtx_mva  = l.vtx_std_evt_mva->at(diphoton_id);
      float sigmaMrv = massResolutionCalculator->massResolutionEonly();
      float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
      float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
      // easy to calculate vertex probability from vtx mva output
      float vtxProb   = 1.-0.49*(vtx_mva+1.0); /// should better use this: vtxAna_.setPairID(diphoton_id); vtxAna_.vertexProbability(vtx_mva); PM

      float phoid_mvaout_lead = l.photonIDMVANew(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4,bdtTrainingPhilosophy.c_str());
      float phoid_mvaout_sublead = l.photonIDMVANew(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4,bdtTrainingPhilosophy.c_str());

      // apply di-photon level smearings and corrections
      selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);

      if( cur_type != 0 && doMCSmearing ) {
	applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, phoid_mvaout_lead,phoid_mvaout_sublead,diPhoSys, syst_shift);
      }
                           
      // Must be calculated after photon id has potentially been smeared
      diphobdt_output = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id] ,vtxProb,lead_p4,sublead_p4 ,sigmaMrv,sigmaMwv,sigmaMeonly ,bdtTrainingPhilosophy.c_str() ,phoid_mvaout_lead,phoid_mvaout_sublead);
      kinematic_bdtout = diphobdt_output;

      mass = Higgs.M();
      float ptHiggs = Higgs.Pt();

      bool isEBEB  = fabs(lead_p4.Eta() < 1.4442 ) && fabs(sublead_p4.Eta()<1.4442);
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
    if( mass>=massMin && mass<=massMax) l.FillHist2D("passVBFCiC_vs_passVBFMVA",0, float(VBFevent_MVA), float(VBFevent_CIC), weight*genLevWeight);

    return false;
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
	    int cat = categoryFromBoundaries( bdtCategoryBoundaries, bdtout );
	    if( VBFevent && cat > -1 ) cat = bdtCategoryBoundaries.size();
	    return cat;
    } else std::cerr << "No BDT Philosophy known - " << bdtTrainingPhilosophy << std::endl;
}

void CategoryAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    if ( doEresolSmear ) eResolSmearer->resetRandom();
}



// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
