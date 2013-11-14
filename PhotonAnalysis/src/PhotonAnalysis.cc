#include "../interface/PhotonAnalysis.h"
#include "../interface/BTagWeight.h"


#include "PhotonReducedInfo.h"
#include "Sorters.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "JetAnalysis/interface/JetHandler.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include "TH2F.h"

#define PADEBUG 0

using namespace std;

// Return the cosine of the angle between the momentums of 2 TLorentzVectors
Double_t CosAngle(TLorentzVector& p4_1, TLorentzVector& p4_2)
{
  Double_t ptot2 = p4_1.Vect().Mag2()*p4_2.Vect().Mag2();
  if(ptot2 <= 0)
    return 0.0;
  else
  {
    Double_t cosTheta = p4_1.Vect().Dot(p4_2.Vect())/TMath::Sqrt(ptot2);
    if(cosTheta > 1.0)
      cosTheta = 1.0;
    if(cosTheta < -1.0)
      cosTheta = -1.0;
    return cosTheta;
  }
}

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

    reRunVtx = false;
    rematchConversions = true;

    keepPP = true;
    keepPF = true;
    keepFF = true;

    selectprocess = false;
    processtoselect = -1;

    doSystematics = true;

    zero_ = 0.;

    recomputeBetas = false;
    recorrectJets = false;
    rerunJetMva = false;
    recomputeJetWp = false;
    applyJer = false;
    jerShift = 0.;
    applyJecUnc = false;
    jecShift = 0.;
    emulateJetResponse = false;
    jetResponseLumiStep = 0.1;
    jetHandler_ = 0;

    shiftBtagEffUp_bc=0;
    shiftBtagEffDown_bc=0;
    shiftBtagEffUp_l=0;
    shiftBtagEffDown_l=0;

    applyBtagSF=0;
    applyLeptonSF=0;
    nBMC=0;
    nCMC=0;
    nLMC=0;
    nBtagB=0;
    nBtagC=0;
    nBtagL=0;

    ptgg_btag_cut=70.;
    costhetastar_btag_cut=0.57;
    ptjet_btag_cut=30.;

    ptgg_0tag_cut=90.;
    costhetastar_0tag_cut=0.57;
    ptjet_0tag_cut=30.;

    deltaRPholep_cut=0.5;
    doMinvCut=false;
    ptJets_ttH_thresh=30.;
    njets_tthHad_thresh=5;

    removeBtagtth=false;
    createCS=false;

    diphobdt_output_Cut_TTHlep=-100;
    diphobdt_output_Cut_TTHhad=-100;
    diphobdt_output_Cut_VHhadBtag=-100;
    diphobdt_output_Cut_VHhad=-100;


    optimizeMVA=false;

    doLooseLep=false;
    doDrGsfTrackCut=false;
    ptjet_loosecut=25.;

    drSC_lep=0.;
    drGsf_lep=0.;


    reComputeCiCPF = false;
    skimOnDiphoN = true;

    doStocasticSmearingSyst = false;
    splitEscaleSyst = false;
    scale_offset_corr_error_file = "";
    splitEresolSyst = false;
    corr_smearing_file = "";
    mass_resol_file = "";

    run7TeV4Xanalysis = false;

    emulateBeamspot = false;
    reweighBeamspot = false;
    saveBSTrees_=false;
    beamspotWidth   = 0.;
    emulatedBeamspotWidth = 0.;
    targetsigma=1.0;
    sourcesigma=1.0;
    rescaleDZforVtxMVA=false;
    pickRandomVtx = false;
    randomizeHiggsPt = false;

    saveDatacardTrees_=false;
    saveSpinTrees_=false;
    saveVBFTrees_=false;
    
    doDiphoMvaUpFront=false;
    useGbrVbfMva=false;
    useGbrDiphotonMva=false;
    
    mvaVbfSelection=false;
    mvaVbfUseDiPhoPt=true;
    mvaVbfUsePhoPt=true;
    mvaVbfSpin=false;
    bookDiPhoCutsInVbf=false;
    multiclassVbfSelection=false;
    vbfVsDiphoVbfSelection=false;
    
    combinedmvaVbfSelection=false;
    reweighPt=false;

    _foresteb=0;  
    _forestee=0; 
}

// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::~PhotonAnalysis()
{
    if( jetHandler_ != 0 ) delete jetHandler_;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Term(LoopAll& l)
{}

std::ostream & operator<< (std::ostream & out, const PhotonCategory & pcat) 
{
    out << pcat.name << ": "  
	<< " (" << pcat.minet  << "," << pcat.maxet  << ") " 
	<< " (" << pcat.mineta << "," << pcat.maxeta << ") "
	<< " (" << pcat.minr9  << "," << pcat.maxr9  << ") "
	<< " " << pcat.type
	;
}

// ----------------------------------------------------------------------------------------------------
void readEnergyScaleOffsets(const std::string &fname, EnergySmearer::energySmearingParameters::eScaleVector &escaleOffsets,
                            EnergySmearer::energySmearingParameters::phoCatVector &photonCategories, bool data=true
                            )
{
    // read in energy scale corrections to be applied in run ranges
    std::cout << "readEnergyScaleOffsets " << fname << std::endl;
    std::fstream in(fname.c_str());
    assert( in );
    char line[200];
    float EBHighR9, EBLowR9, EBm4HighR9, EBm4LowR9, EEHighR9, EELowR9;
    char catname[200];
    do {
        int nread = 0;
        float minet=0., maxet=1.e+9, mineta=0., maxeta=999., minr9=-999, maxr9=999, offset=0., stocastic=0., err=0., stocastic_err=0., 
            pivot=0.;
        int type;
        int  first, last;
        
        in.getline( line, 200, '\n' );
        if( line[0] == '#' ) { continue; }
        
        if( (nread == 0) && 
            (nread = sscanf(line,"%s %d %f %f %f %f %f %f %d %d %f %f %f %f %f\n", &catname, &type, 
                            &minet, &maxet, &mineta, &maxeta, &minr9, &maxr9, &first, &last, &pivot, &offset, &err, 
                            &stocastic, &stocastic_err ) ) != 15 ) {
            nread=0, minet=0., maxet=1.e+9, mineta=0., maxeta=999., minr9=-999, maxr9=999, offset=0., stocastic=0., stocastic_err=0.,
                err=0., pivot=0;
        }
        if( (nread == 0) && 
            (nread = sscanf(line,"%s %d %f %f %f %f %d %d %f %f %f %f %f\n", &catname, &type, 
                            &mineta, &maxeta, &minr9, &maxr9, &first, &last, &pivot, &offset, &err, 
                            &stocastic, &stocastic_err ) ) != 13 ) {
            nread=0, minet=0., maxet=1.e+9, mineta=0., maxeta=999., minr9=-999, maxr9=999, offset=0., stocastic=0., stocastic_err=0.,
                err=0., pivot=0;
        }
        if( (nread == 0) &&
            ( nread = sscanf(line,"%s %d %f %f %f %f %d %d %f %f\n", &catname, &type, 
                             &mineta, &maxeta, &minr9, &maxr9, &first, &last, &offset, &err  ) ) != 10 ) {
            nread=0, minet=0., maxet=1.e+9, mineta=0., maxeta=999., minr9=-999, maxr9=999, offset=0., stocastic=0., stocastic_err=0.,
                err=0., pivot=0;
        }
        if( nread == 0 ) { 
            continue; 
        }
	
        std::cout << "Energy scale (or smearing) by run " << nread  << " " << catname   << " " << type 
                  << " " << minet << "<E_T<" << maxet << " " << mineta << "<eta<" << maxeta    << " " << minr9 << "<r9<" << maxr9 
                  << " " << first << "<run<" << last  << " " << pivot << " " << offset << "+-" << err 
                  << " " << stocastic << "+-" << stocastic_err << std::endl;
        
        assert( type>=0 && type<=2 );
        
        EnergySmearer::energySmearingParameters::eScaleVector::reverse_iterator escaleOffset =
            find(escaleOffsets.rbegin(),escaleOffsets.rend(),std::make_pair(first,last));
        if( escaleOffset == escaleOffsets.rend() ) {
            std::cout << "  adding new range range " << first << " " << last << std::endl;
            escaleOffsets.push_back(EnergyScaleOffset(first,last));
            escaleOffset = escaleOffsets.rbegin();
        }
        // chck if the category is already defined
        if( find(photonCategories.begin(), photonCategories.end(), std::string(catname) ) == photonCategories.end() ) {
            photonCategories.push_back(PhotonCategory(minet,maxet,mineta,maxeta,minr9,maxr9,(PhotonCategory::photon_type_t)type,catname));
            std::cout << "  defining new category " << photonCategories.back() << std::endl;
        }
        // assign the scale offset and error for this category and this run range
        escaleOffset->scale_offset[catname] = data ? -offset : offset;
        escaleOffset->scale_stocastic_offset[catname] = data ? -stocastic : stocastic;
        escaleOffset->scale_offset_error[catname] = err;
        escaleOffset->scale_stocastic_offset_error[catname] = stocastic_err;
        escaleOffset->scale_stocastic_pivot[catname] = pivot;
        
    } while( in );
    
    in.close();
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuMap(const char * fname, TDirectory * dir, TH1 * target)
{
    std::fstream in(fname);
    assert( in );
    char line[200];
    int typid;
    char dname[200];
    do {
        in.getline( line, 200, '\n' );

        if( sscanf(line,"%d %s",&typid,dname) != 2 ) { continue; }
        std::cout << "Reading PU weights for sample " << typid << " from " << dname << std::endl;
        TDirectory * subdir = (TDirectory *)dir->Get(dname);
        assert( subdir != 0 );
        cout << "loadPuMap " << typid << " from subdir " << endl;
        loadPuWeights(typid, subdir, target);
    } while ( in );
    in.close();
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuWeights(int typid, TDirectory * puFile, TH1 * target)
{
    cout<<"Loading pileup weights for typid " << typid <<endl;

    TH1 * hweigh = (TH1*) puFile->Get("pileup_weights");
    if( hweigh == 0 ) {
	  if( PADEBUG ) std::cout << "pileup_weights not found: trying weights " << endl;
	  hweigh = (TH1*) puFile->Get("weights");
    }

    TH1 * gen_pu = (TH1*)puFile->Get("pileup");
    if( gen_pu == 0 ) {
	if( PADEBUG ) std::cout << "pileup not found: trying generated_pu " << endl;
	gen_pu = (TH1*)puFile->Get("generated_pu");
    }

    assert( gen_pu != 0 );
    bool delHwei = false;
    if( target != 0 ) {
	// compute weights on the fly based on the target pileup scenario
	if( hweigh == 0 ) {
	    hweigh = (TH1*)gen_pu->Clone();
	    delHwei = true;
	}
	hweigh->Reset("ICE");
	for( int ii=1; ii<hweigh->GetNbinsX(); ++ii ) {
	    hweigh->SetBinContent( ii, target->GetBinContent( target->FindBin( hweigh->GetBinCenter(ii) ) ) );
	}
	hweigh->Divide(hweigh, gen_pu, 1., 1./gen_pu->Integral() );
    } else {
	// Normalize weights such that the total cross section is unchanged
	TH1 * eff = (TH1*)hweigh->Clone("eff");
	eff->Multiply(gen_pu);
	hweigh->Scale( gen_pu->Integral() / eff->Integral()  );
	delete eff;
    }
    
    weights[typid].clear();
    for( int ii=1; ii<hweigh->GetNbinsX(); ++ii ) {
	weights[typid].push_back(hweigh->GetBinContent(ii));
    }

    std::cout << "pile-up weights: ["<<typid<<"]";
    std::copy(weights[typid].begin(), weights[typid].end(), std::ostream_iterator<double>(std::cout,","));
    std::cout << std::endl;

    if( delHwei ) { delete hweigh; }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::load2DPuWeights(int typid, TDirectory * puFile, std::vector<TH1*> targets) {

    cout<<"Loading 2D pileup weights for typid " << typid << " for all run periods." << endl;

    rd_weights[typid].resize(targets.size());    
    for (unsigned int i=0; i<targets.size(); i++) {
	std::cout<<targets[i]->GetName() << std::endl;
	TH1* gen_pu = ((TH2F*)puFile->Get("pu_2D"))->ProjectionX("pileup", i+1, i+1);
	assert( gen_pu != 0 );
	TH1* hweigh =(TH1*)targets[i]->Clone();
	hweigh->Reset("ICE");
	hweigh->Divide(targets[i], gen_pu, 1., 1./gen_pu->Integral());
	rd_weights[typid][i].clear();
	for( int ii=1; ii<hweigh->GetNbinsX(); ++ii ) {
	    rd_weights[typid][i].push_back(hweigh->GetBinContent(ii));
	}
    }

    std::cout << "pile-up 2D weights: ["<<typid<<"]";
    //std::cout << targets.size() << std::endl;
    for (int i=0; i<targets.size(); i++) {
	std::copy(rd_weights[typid][i].begin(),rd_weights[typid][i].end(), std::ostream_iterator<double>(std::cout,","));
	std::cout << std::endl;
    }
}

// ----------------------------------------------------------------------------------------------------
float PhotonAnalysis::getPuWeight(int n_pu, int sample_type, SampleContainer* container, bool warnMe, int run) {

    if ( sample_type !=0 && puHist != "") {

        bool hasSpecificWeight = weights.find( sample_type ) != weights.end();
	bool hasSpecificWeight2D = rd_weights.find(sample_type) != rd_weights.end();
	
	if(!hasSpecificWeight && !hasSpecificWeight2D && container != 0 && container->pileup != 0 ) {
	    std::cout << "On-the fly pileup reweighing typeid " << sample_type << " " << container->pileup << std::endl;
	    TFile * samplePu = TFile::Open(container->pileup.c_str());
	    TKey* key = samplePu->FindKey("pu_2D");
	    if (key != 0 && !puTargets.empty()) {
		load2DPuWeights(sample_type, samplePu, puTargetHists);
		hasSpecificWeight2D = true;
	    } else {
		loadPuWeights(sample_type, samplePu, puTargetHist);
		hasSpecificWeight = true;
	    }
	    samplePu->Close();
	} else if( sample_type < 0 && !hasSpecificWeight && !hasSpecificWeight2D && warnMe) {
            std::cerr  << "WARNING no pu weights specific for sample " << sample_type << std::endl;
        }

	if (hasSpecificWeight2D) {
	    Int_t index = -1;
	    
	    if (run <= 197495)
		index = 0;
	    else if (run > 197495 && run <= 203767)
		index = 1;
	    else if (run > 203767)
		index = 2;
	    else {
		std::cout << "ERROR: it is not possible to weight this event (unrecognized run " << run << ")" << std::endl;
		abort();
	    }
	    std::vector<double>& puweights = rd_weights[sample_type][index];
	    if (n_pu < puweights.size()) {
		return puweights[n_pu];
	    }
	} 

	if (hasSpecificWeight) {
	    std::vector<double> & puweights = weights[sample_type];
	    if(n_pu<puweights.size()){
		return puweights[n_pu];
	    }
	    else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
		cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< sample_type <<"], event will not be reweighted for pileup"<<endl;
	    }
	}
    }

    return 1.;
}

// ----------------------------------------------------------------------------------------------------
//float PhotonAnalysis::getPuWeight(int n_pu, int sample_type, SampleContainer* container, bool warnMe)
//{
//    if ( sample_type !=0 && puHist != "") {
//        bool hasSpecificWeight = weights.find( sample_type ) != weights.end() ;
//    if( ! hasSpecificWeight && container != 0 && container->pileup != 0 ) {
//        std::cout << "On-the fly pileup reweighing typeid " << sample_type << " " << container->pileup << std::endl;
//        TFile * samplePu = TFile::Open(container->pileup.c_str());
//        loadPuWeights(sample_type, samplePu, puTargetHist);
//        samplePu->Close();
//        hasSpecificWeight = true;
//    } else if( sample_type < 0 && !hasSpecificWeight && warnMe ) {
//            std::cout  << "WARNING no pu weights specific for sample " << sample_type << std::endl;
//        }
//        std::vector<double> & puweights = hasSpecificWeight ? weights[ sample_type ] : weights[0];
//        if(n_pu<puweights.size()){
//            return puweights[n_pu];
//        }
//        else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
//            cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< sample_type <<"], event will not be reweighted for pileup"<<endl;
//        }
//    }
//    return 1.;
//}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::applyGenLevelSmearings(double & genLevWeight, const TLorentzVector & gP4, int npu, int sample_type, BaseGenLevelSmearer * sys, float syst_shift)
{
    static int nwarnings=10;
    for(std::vector<BaseGenLevelSmearer*>::iterator si=genLevelSmearers_.begin(); si!=genLevelSmearers_.end(); si++){
    float genWeight=1;
    if( sys != 0 && *si == *sys ) {
        (*si)->smearEvent(genWeight, gP4, npu, sample_type, syst_shift );
    } else {
        (*si)->smearEvent(genWeight, gP4, npu, sample_type, 0. );
    }
    if( genWeight < 0. ) {
        if( syst_shift == 0. ) {
        std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
        assert(0);
        } else {
        if( nwarnings-- > 0 ) {
            std::cout <<  "WARNING: negative during systematic scan in " << (*si)->name() << " " << genWeight << std::endl;
        }
        genWeight = 0.;
        }
    }
    genLevWeight*=genWeight;
    }
    // Set some provate members with info
    generatorPt_ = gP4.Pt();
    generatorY_  = gP4.Rapidity();
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::applySinglePhotonSmearings(std::vector<float> & smeared_pho_energy, std::vector<float> & smeared_pho_r9, std::vector<float> & smeared_pho_weight,
                                                int cur_type, const LoopAll & l, const float * energyCorrected, const float * energyCorrectedError,
                                                BaseSmearer * sys, float syst_shift

    )
{
    static int nwarnings = 10;
    static int cache_run = -1, cache_lumis = -1, cache_event = -1;
    bool fillInfo = false;
    if( l.run != cache_run || l.lumis != cache_lumis || l.event != cache_event ) {
        fillInfo = true;
        cache_run   = l.run  ;
        cache_lumis = l.lumis;
        cache_event = l.event;
    }

    /// if( sys ) std::cout << "applySinglePhotonSmearings " << fillInfo << " " << syst_shift << " " << ( sys != 0 ? sys->name() : " " ) <<std::endl;

    if( fillInfo ) {
        photonInfoCollection.clear();
    }
    smeared_pho_energy.resize(l.pho_n,0.);
    smeared_pho_r9.resize(l.pho_n,0.);
    smeared_pho_weight.resize(l.pho_n,0.);
    for(int ipho=0; ipho<l.pho_n; ++ipho ) {
        
        std::vector<std::vector<bool> > p;
        if( fillInfo ) {
            PhotonReducedInfo info (
                *((TVector3*)     l.sc_xyz->At(l.pho_scind[ipho])),
                ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(),
                energyCorrected[ipho],
                l.pho_isEB[ipho], l.pho_r9[ipho],                
                ipho,
                true, // WARNING  setting pass photon ID flag for all photons. This is safe as long as only selected photons are used
                (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
                );
            photonInfoCollection.push_back(info);
        } else {
            photonInfoCollection[ipho].reset();
        }
        PhotonReducedInfo & phoInfo = photonInfoCollection[ipho];
        //// if( sys ) phoInfo.dump();
        
        int ieta, iphi;
        l.getIetaIPhi(ipho,ieta,iphi);
        phoInfo.addSmearingSeed( (unsigned int)l.sc_raw[l.pho_scind[ipho]] + abs(ieta) + abs(iphi) + l.run + l.event + l.lumis );
        phoInfo.setSphericalPhoton(l.CheckSphericalPhoton(ieta,iphi));
        
        // FIXME add seed to syst smearings
        
        float pweight = 1.;
        // smear MC. But apply energy corrections and scale adjustement to data
        if( cur_type != 0 && doMCSmearing ) {
            for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
                float sweight = 1.;
                if( sys != 0 && *si == *sys ) {
                    // move the smearer under study by syst_shift
                    (*si)->smearPhoton(phoInfo,sweight,l.run,syst_shift);
                } else {
                    // for the other use the nominal points
                    (*si)->smearPhoton(phoInfo,sweight,l.run,0.);
                }
                if( sweight < 0. ) {
                    if( syst_shift == 0. ) {
                        std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
                        assert(0);
                    } else {
                        if( nwarnings-- > 0 ) {
                            std::cout <<  "WARNING: negative during systematic scan in " << (*si)->name() << std::endl;
                        }
                        sweight = 0.;
                    }
                }
                pweight *= sweight;
            }
        } else if( cur_type == 0 ) {
            float sweight = 1.;
            if( doEcorrectionSmear )  {
                eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
            }
            eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
            pweight *= sweight;
        }
        
        //// if( ipho == 0 ) { phoInfo.dump(); }
        smeared_pho_energy[ipho] = phoInfo.energy();
        smeared_pho_r9[ipho]     = phoInfo.r9();
        smeared_pho_weight[ipho] = pweight;

    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::fillDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs,
                  float & lead_r9, float & sublead_r9, TVector3 *& vtx, const float * energy,
                  const LoopAll & l, int diphoton_id)
{
    fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, energy,
                  l, l.dipho_leadind[diphoton_id], l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id]);
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::fillDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs,
                  float & lead_r9, float & sublead_r9, TVector3 *& vtx, const float * energy,
                  const LoopAll & l, int leadind, int subleadind, int myvtx)
{
    int vtx_ind = myvtx;

    int cur_type = l.itype[l.current];
    bool isCorrectVertex;
    if (cur_type!=0) isCorrectVertex= (*((TVector3*)l.vtx_std_xyz->At(vtx_ind))-*((TVector3*)l.gv_pos->At(0))).Mag() < 1.;

    TLorentzVector lead_p4_bef = l.get_pho_p4( leadind, vtx_ind, energy);
    TLorentzVector sublead_p4_bef = l.get_pho_p4( subleadind, vtx_ind, energy);
    TLorentzVector Higgs_bef = lead_p4_bef + sublead_p4_bef;

    //TRandom3 rand;
    //rand.SetSeed(0);
    bool changed=false;
    if (emulateBeamspot && cur_type!=0 && (lastRun!=l.run || lastEvent!=l.event || lastLumi!=l.lumis)){
        double genVtxZ = ((TVector3*)l.gv_pos->At(0))->Z();
        double myVtxZ = ((TVector3*)l.vtx_std_xyz->At(vtx_ind))->Z();
        ((TVector3*)l.vtx_std_xyz->At(vtx_ind))->SetZ(genVtxZ+(targetsigma/sourcesigma)*(myVtxZ-genVtxZ));
        changed=true;
    }
    
    if (pickRandomVtx && cur_type!=0 && (lastRun!=l.run || lastEvent!=l.event || lastLumi!=l.lumis)){
        TRandom3 rand;
        rand.SetSeed(0);
        double randVtxZ = rand.Gaus(0.,4.8);
        ((TVector3*)l.vtx_std_xyz->At(vtx_ind))->SetZ(randVtxZ);
    }

    lead_p4 = l.get_pho_p4( leadind, vtx_ind, energy);
    sublead_p4 = l.get_pho_p4( subleadind, vtx_ind, energy);
    lead_r9    = l.pho_r9[leadind];
    sublead_r9 = l.pho_r9[subleadind];
    Higgs = lead_p4 + sublead_p4;
    vtx = (TVector3*)l.vtx_std_xyz->At(vtx_ind);
    
    if (randomizeHiggsPt && cur_type!=0 && (lastRun!=l.run || lastEvent!=l.event || lastLumi!=l.lumis)){
        TRandom3 rand;
        rand.SetSeed(0);
        Higgs.SetPerp(Higgs.Pt()*rand.Gaus(1.,0.1));
    }

    if (changed && PADEBUG){
        cout << "emBS: " << emulateBeamspot << " rwBS: " << reweighBeamspot << " tS: " << targetsigma << " sS: " << sourcesigma << " bW: " << beamspotWidth << endl;
        cout << "Before (eta1,eta2,mh): " << lead_p4_bef.Eta() << " " << sublead_p4_bef.Eta() << " " << Higgs_bef.M() << endl;
        cout << "After  (eta1,eta2,mh): " << lead_p4.Eta() << " " << sublead_p4.Eta() << " " << Higgs.M() << endl;
    }

    lastRun = l.run;
    lastEvent = l.event;
    lastLumi = l.lumis;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::applyDiPhotonSmearings(TLorentzVector & Higgs, TVector3 & vtx, int category, int cur_type, const TVector3 & truevtx,
                        float & evweight, float & idmva1, float & idmva2,
                        BaseDiPhotonSmearer * sys, float syst_shift)
{
    static int nwarnings=10;
    float pth = Higgs.Pt();
    for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
        float rewei=1.;
	if( sys != 0 && *si == *sys ) {
	    (*si)->smearDiPhoton( Higgs, vtx, rewei, category, cur_type, truevtx, idmva1, idmva2, syst_shift );
	} else {
	    (*si)->smearDiPhoton( Higgs, vtx, rewei, category, cur_type, truevtx, idmva1, idmva2, 0. );
	}
	if( rewei < 0. ) {
	    if( syst_shift == 0. ) {
		std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		assert(0);
	    } else {
		if( nwarnings-- > 0 ) {
		    std::cout <<  "WARNING: negative during systematic scan in " << (*si)->name() << std::endl;
		}
		rewei = 0.;
	    }
	}
	evweight *= rewei;
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Init(LoopAll& l)
{
    if(PADEBUG)
        cout << "InitRealPhotonAnalysis START"<<endl;

    if(energyCorrectionMethod=="DaunceyAndKenzie"){
        energyCorrected     = (l.pho_residCorrEnergy);
        energyCorrectedError= (l.pho_residCorrResn);
    }else if(energyCorrectionMethod=="Bendavid"){
        energyCorrected     = (l.pho_regr_energy);
        energyCorrectedError= (l.pho_regr_energyerr);
    }else if(energyCorrectionMethod=="BendavidOTF"){
        energyCorrected     = (l.pho_regr_energy_otf);
        energyCorrectedError= (l.pho_regr_energyerr_otf);
    }else{
        assert(doEcorrectionSmear==false);
    }
    if (doEcorrectionSmear) std::cout << "using energy correction type: " << energyCorrectionMethod << std::endl;
    else                    std::cout << "NOT using energy correction (sbattogiu)"<< std::endl;

    if( vtxVarNames.empty() ) {
        vtxVarNames.push_back("ptbal"), vtxVarNames.push_back("ptasym"), vtxVarNames.push_back("logsumpt2");
    }

    /// // trigger

    if (l.runZeeValidation && l.sqrtS != 7) {
        
        triggerSelections.push_back(TriggerSelection(1,-1));
        if (useRunDTriggersForZee) {
            //use this if the validation includes runD
            triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v");
        } else {
            //use this if the validation does not include runD
            triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v");
        }
        
        triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
        
    } else {
        
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
        
        triggerSelections.push_back(TriggerSelection(194270,203767));
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
        
        triggerSelections.push_back(TriggerSelection(203768,-1));
        triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v");
        triggerSelections.back().addpath("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v");
        triggerSelections.back().addpath("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
        triggerSelections.back().addpath("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
        triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
        triggerSelections.back().addpath("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v");
        triggerSelections.back().addpath("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
        triggerSelections.back().addpath("HLT_Photon36_R9Id85_Photon22_R9Id85_v");
    }
    
    if( l.typerun != l.kReduce ) {
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
        
        if( mvaVbfSelection || multiclassVbfSelection || combinedmvaVbfSelection ) {
          
          if( useGbrVbfMva ) {
              /// jeteta1,jeteta2,jetpt1,jetpt2,zeppenfeld,dphidijetgg,dijetmass,ptgg
              TFile * fin = TFile::Open(gbrVbfFile.c_str());
              RooWorkspace * ws = (RooWorkspace*) (fin->Get("wsfitmc")->Clone());
              
              gbrVbfReader_ = new RooFuncReader(ws,"sigVBFxxb","trainingvars");
              gbrVbfReader_->bookVariable("jeteta1",    &myVBFLeadJEta);
              gbrVbfReader_->bookVariable("jeteta2",    &myVBFSubJEta);
              gbrVbfReader_->bookVariable("jetpt1",     &myVBFLeadJPt);
              gbrVbfReader_->bookVariable("jetpt2",     &myVBFSubJPt);
              gbrVbfReader_->bookVariable("zeppenfeld", &myVBFZep);
              gbrVbfReader_->bookVariable("dphidijetgg",&myVBFdPhiTrunc);
              gbrVbfReader_->bookVariable("dijetmass",  &myVBF_Mjj);	   
              gbrVbfReader_->bookVariable("ptgg",       &myVBFDiPhoPtOverM);
              fin->Close();

            if( combinedmvaVbfSelection ) {
                fin = TFile::Open(gbrVbfDiPhoFile.c_str());
                ws = (RooWorkspace*) (fin->Get("wsfitmc")->Clone());
                gbrVbfDiphoReader_ = new RooFuncReader(ws,"sigVBFxxb","trainingvars");
                gbrVbfDiphoReader_->bookVariable("masserr",         l.tmva_dipho_MIT_dmom);
                gbrVbfDiphoReader_->bookVariable("masserrwrong",    l.tmva_dipho_MIT_dmom_wrong_vtx);
                gbrVbfDiphoReader_->bookVariable("vtxprob",         l.tmva_dipho_MIT_vtxprob);
                gbrVbfDiphoReader_->bookVariable("pt1",             l.tmva_dipho_MIT_ptom1);
                gbrVbfDiphoReader_->bookVariable("pt2",             l.tmva_dipho_MIT_ptom2);
                gbrVbfDiphoReader_->bookVariable("eta1",            l.tmva_dipho_MIT_eta1);
                gbrVbfDiphoReader_->bookVariable("eta2",            l.tmva_dipho_MIT_eta2);
                gbrVbfDiphoReader_->bookVariable("dphi",            l.tmva_dipho_MIT_dphi);
                gbrVbfDiphoReader_->bookVariable("idmva1",          l.tmva_dipho_MIT_ph1mva);
                gbrVbfDiphoReader_->bookVariable("idmva2",          l.tmva_dipho_MIT_ph2mva);
                gbrVbfDiphoReader_->bookVariable("jeteta1",         &myVBFLeadJEta);
                gbrVbfDiphoReader_->bookVariable("jeteta2",         &myVBFSubJEta);
                gbrVbfDiphoReader_->bookVariable("jetpt1",          &myVBFLeadJPt);
                gbrVbfDiphoReader_->bookVariable("jetpt2",          &myVBFSubJPt);
                gbrVbfDiphoReader_->bookVariable("zeppenfeld",      &myVBFZep);
                gbrVbfDiphoReader_->bookVariable("dphidijetgg",     &myVBFdPhiTrunc);
                gbrVbfDiphoReader_->bookVariable("dijetmass",       &myVBF_Mjj);	   
                fin->Close();
            }
            
          } else {
            tmvaVbfReader_ = new TMVA::Reader( "!Color:!Silent" );
            
            if (combinedmvaVbfSelection) {
              tmvaVbfReader_->AddVariable("dijet_leadEta",    &myVBFLeadJEta);
              tmvaVbfReader_->AddVariable("dijet_subleadEta", &myVBFSubJEta);
              tmvaVbfReader_->AddVariable("dijet_LeadJPt",    &myVBFLeadJPt);
              tmvaVbfReader_->AddVariable("dijet_SubJPt",     &myVBFSubJPt);
              tmvaVbfReader_->AddVariable("dijet_Zep",        &myVBFZep);
              if (!combinedmvaVbfSelection) {
                tmvaVbfReader_->AddVariable("dijet_dPhi",       &myVBFdPhi);
              } else {
                if(l.sqrtS==7){
                  tmvaVbfReader_->AddVariable("min(dijet_dPhi,2.9416)", &myVBFdPhiTrunc);
                } else if(l.sqrtS==8){
                  tmvaVbfReader_->AddVariable("min(dijet_dPhi,2.916)", &myVBFdPhiTrunc);
                } else {
                  std::cout<<"sqrtS is not 7 or 8 but is "<<l.sqrtS<<std::endl;
                }
              }
              tmvaVbfReader_->AddVariable("dijet_Mjj",        &myVBF_Mjj);	   
              tmvaVbfReader_->AddVariable("dipho_pt/mass",    &myVBFDiPhoPtOverM);
              
              tmvaVbfDiphoReader_ = new TMVA::Reader("!Color:!Silent"); 
              //tmvaVbfDiphoReader_->AddVariable("bdt_incl",                       &myVBFDIPHObdt);
              //tmvaVbfDiphoReader_->AddVariable("bdt_dijet_sherpa_plusdiphoptom", &myVBF_MVA);
              //tmvaVbfDiphoReader_->AddVariable("dipho_pt/mass",                  &myVBFDiPhoPtOverM);
              tmvaVbfDiphoReader_->AddVariable("dipho_mva",                       &myVBFDIPHObdt);
              if(l.sqrtS==7){
                tmvaVbfDiphoReader_->AddVariable("bdt_dijet_7TeV_ptrewght",         &myVBF_MVA);
              } else if(l.sqrtS==8){
                tmvaVbfDiphoReader_->AddVariable("bdt_dijet_maxdPhi",               &myVBF_MVA);
              } else {
                std::cout<<"sqrtS is not 7 or 8 but is "<<l.sqrtS<<std::endl;
              }
              tmvaVbfDiphoReader_->AddVariable("dipho_pt/mass",                  &myVBFDiPhoPtOverM);
              tmvaVbfDiphoReader_->BookMVA(mvaVbfDiphoMethod, mvaVbfDiphoWeights);
            } else {
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
            }
            
            tmvaVbfReader_->BookMVA( mvaVbfMethod, mvaVbfWeights );
          }
        }
        
        if( mvaVbfSpin && (mvaVbfSelection || multiclassVbfSelection) )
        {
            tmvaVbfSpinReader_ = new TMVA::Reader( "!Color:!Silent" );
            
            tmvaVbfSpinReader_->AddVariable("absDeltaPhiJJ := abs(deltaPhiJJ)", &myVBFSpin_absDeltaPhiJJ);
            tmvaVbfSpinReader_->AddVariable("absCosThetaJ1 := abs(cosThetaJ1)", &myVBFSpin_absCosThetaJ1);
            tmvaVbfSpinReader_->AddVariable("absCosThetaJ2 := abs(cosThetaJ2)", &myVBFSpin_absCosThetaJ2);
            
            //tmvaVbfSpinReader_->AddVariable("absDeltaPhiJJS := abs(deltaPhiJJS)", &myVBFSpin_absDeltaPhiJJS);
            tmvaVbfSpinReader_->AddVariable("absCosThetaS := abs(cosThetaS)", &myVBFSpin_absCosThetaS);
            tmvaVbfSpinReader_->AddVariable("absDeltaPhiJJL := abs(deltaPhiJJL)", &myVBFSpin_absDeltaPhiJJL);
            tmvaVbfSpinReader_->AddVariable("absCosThetaL := abs(cosThetaL)", &myVBFSpin_absCosThetaL);
            
            tmvaVbfSpinReader_->BookMVA( mvaVbfSpinMethod, mvaVbfSpinWeights );
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
    readEnergyScaleOffsets(scale_offset_file, eSmearDataPars.scale_offset_byrun, eSmearDataPars.photon_categories);
    // if the scale offset file defines the categories set the category type to automatic
    if( ! eSmearDataPars.photon_categories.empty() ) {
        eSmearDataPars.categoryType = "Automagic";
        eSmearDataPars.n_categories = -1;
    }

    // energy scale corrections to Data
    eScaleDataSmearer = new EnergySmearer( eSmearDataPars );
    eScaleDataSmearer->name("E_scale_data");
    eScaleDataSmearer->doEnergy(true);
    eScaleDataSmearer->scaleOrSmear(true);

    // Read energy scale errors and energy smaerings from dat files
    assert( ! scale_offset_error_file.empty() && ! smearing_file.empty() );
    
    // Use the same format used for the run-dependent energy corrections
    EnergySmearer::energySmearingParameters::eScaleVector tmp_scale_offset, tmp_smearing;
    EnergySmearer::energySmearingParameters::phoCatVector tmp_scale_cat, tmp_smearing_cat;
    readEnergyScaleOffsets(scale_offset_error_file, tmp_scale_offset, tmp_scale_cat,false);
    readEnergyScaleOffsets(smearing_file, tmp_smearing, tmp_smearing_cat,false);
    
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
    eSmearPars.scale_stocastic_offset = tmp_scale_offset[0].scale_stocastic_offset;
    eSmearPars.scale_stocastic_offset_error = tmp_scale_offset[0].scale_stocastic_offset_error;
    eSmearPars.scale_stocastic_pivot = tmp_scale_offset[0].scale_stocastic_pivot;
    
    eSmearPars.smearing_sigma = tmp_smearing[0].scale_offset;
    eSmearPars.smearing_sigma_error = tmp_smearing[0].scale_offset_error;
    eSmearPars.smearing_stocastic_sigma = tmp_smearing[0].scale_stocastic_offset;
    eSmearPars.smearing_stocastic_sigma_error = tmp_smearing[0].scale_stocastic_offset_error;
    eSmearPars.smearing_stocastic_pivot = tmp_smearing[0].scale_stocastic_pivot;

    // Energy resolution parameters used for diphotonBDT input
    massResoPars = eSmearPars;
    if( ! mass_resol_file.empty() ) {
        EnergySmearer::energySmearingParameters::eScaleVector tmp_smearing;
        EnergySmearer::energySmearingParameters::phoCatVector tmp_smearing_cat;
        readEnergyScaleOffsets(mass_resol_file, tmp_smearing, tmp_smearing_cat,false);

        // make sure that the scale correction and smearing info is as expected
        assert( tmp_smearing.size() == 1 );
        assert( ! tmp_smearing_cat.empty() );

        // copy the read info to the smarer parameters
        massResoPars.categoryType = "Automagic";
        massResoPars.byRun = false;
        massResoPars.n_categories = tmp_smearing_cat.size();
        massResoPars.photon_categories = tmp_smearing_cat;

        massResoPars.smearing_sigma = tmp_smearing[0].scale_offset;
        massResoPars.smearing_stocastic_sigma = tmp_smearing[0].scale_stocastic_offset;
        massResoPars.smearing_sigma_error = tmp_smearing[0].scale_offset_error;
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

    // MassResolution
    massResolutionCalculator = new MassResolution();

    /* -------------------------------------------------------------------------------------------
       Pileup Reweighting
       https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
       ----------------------------------------------------------------------------------------------  */
        if (puTargets.size() != 0) {
	if (puHist == "auto") {
	    for (unsigned int i=0; i<puTargets.size(); i++) {
		TFile* puTargetFile = TFile::Open(puTargets[i]);
		assert(puTargetFile != 0);
		puTargetHists.push_back((TH1*)puTargetFile->Get("pileup"));
		puTargetHists[i]->SetDirectory(0);
		puTargetHists[i]->Scale(1. / puTargetHists[i]->Integral());
		puTargetFile->Close();
	    }
	} else {
	    std::cout << "WARNING: no other reweighting method implemented for RD MC." << std::endl;
	    abort();
	}
    }

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

    // Jet handling
    if( recomputeBetas || recorrectJets || rerunJetMva || recomputeJetWp || applyJer || applyJecUnc || emulateJetResponse 
	|| l.typerun != l.kFill ) {
	std::cout << "JetHandler: \n"
		  << "recomputeBetas " << recomputeBetas << "\n"
		  << "recorrectJets " << recorrectJets << "\n"
		  << "rerunJetMva " << rerunJetMva << "\n"
		  << "recomputeJetWp " << recomputeJetWp << "\n"
		  << "emulateJetResponse " << emulateJetResponse
		  << std::endl;
	jetHandler_ = new JetHandler(jetHandlerCfg, l);
	jetHandler_->setJetResponseStep(jetResponseLumiStep);
    }

    // Beam spot reweighting
    if( emulateBeamspot || reweighBeamspot ) {
	assert( emulatedBeamspotWidth != 0. );
	beamspotWidth = emulatedBeamspotWidth;
    }
    if( beamspotWidth == 0. ) {
        beamspotWidth = (run7TeV4Xanalysis ? 5.8 : 4.8);
    }
    if (reweighPt) {
        ptreweighfile = TFile::Open(ptreweighfilename.c_str());
        ptreweighHistSM = (TH1F*)ptreweighfile->Get(Form("sm%sRat",ptreweightype.c_str()));
        ptreweighHistGG = (TH1F*)ptreweighfile->Get(Form("gg%sRat",ptreweightype.c_str()));
        ptreweighHistQQ = (TH1F*)ptreweighfile->Get(Form("qq%sRat",ptreweightype.c_str()));
        ptreweighHistSM->SetDirectory(0);
        ptreweighHistGG->SetDirectory(0);
        ptreweighHistQQ->SetDirectory(0);
        ptreweighfile->Close();
    }

    if( l.typerun == LoopAll::kReduce ) {
        // ---------------------- LOAD Regression Classes ---------------------//
        // Implementation copied over from ... 
        // https://github.com/bendavid/GBRLikelihoodEGTools/commit/7aff712aa93c69e5e04664d7556a6bd646af479c#diff-57e3515cb45eaf6857c6bf3a0481aca0
        if (regressionVersion==5){
                //initialize eval vector
            _vals.resize(37);
           
                //load forests from file
                TFile *fgbr = TFile::Open(regressionFile.c_str(),"READ");    
                fgbr->GetObject("EGRegressionForest_EB", _foresteb);
                fgbr->GetObject("EGRegressionForest_EE", _forestee);
                fgbr->Close();

            //recreate pdf with constraint transformations (can't load directly from file due to weird RooWorkspace IO features)
            
            _tgt = new RooRealVar("tgt","",1.);
            _mean = new RooRealVar("mean","",1.);
            _sigma = new RooRealVar("sigma","",1.);
            _n1 = new RooRealVar("n1","",2.);
            _n2 = new RooRealVar("n2","",2.);
            
            _sigmalim = new RooRealConstraint("sigmalim","",*_sigma,0.0002,0.5);
            _meanlim = new RooRealConstraint("meanlim","",*_mean,0.2,2.0);
            _n1lim = new RooRealConstraint("n1lim","",*_n1,1.01,110.);
            _n2lim = new RooRealConstraint("n2lim","",*_n2,1.01,110.);     
            
            _pdf = new RooDoubleCBFast("sigpdf","",*_tgt,RooFit::RooConst(1.),
                           *_sigmalim,RooFit::RooConst(2.0),*_n1lim,RooFit::RooConst(1.0),*_n2lim);
            
            //add to RooArgList for proper garbage collection
            _args.addOwned(*_tgt);
            _args.addOwned(*_mean);
            _args.addOwned(*_sigma);
            _args.addOwned(*_n1);
            _args.addOwned(*_n2);
            _args.addOwned(*_sigmalim);
            _args.addOwned(*_meanlim);
            _args.addOwned(*_n1lim);
            _args.addOwned(*_n2lim);
            _args.addOwned(*_pdf);    
       } 
       else if (regressionVersion==8){ // This is for 7 TeV (we would use V8)
          //initialize eval vector
          _vals.resize(37);
          
          //load forests from file
          if( l.typerun == LoopAll::kReduce ) {
                //load forests from file
                TFile *fgbr = TFile::Open(regressionFile.c_str(),"READ");    
                fgbr->GetObject("EGRegressionForest_EB", _forestDeb);
                fgbr->GetObject("EGRegressionForest_EE", _forestDee);
                fgbr->Close();      
          }
          
          //recreate pdf with constraint transformations (can't load directly from file due to weird RooWorkspace IO features)
          
          _tgt = new RooRealVar("tgt","",1.);
          _mean = new RooRealVar("mean","",1.);
          _sigma = new RooRealVar("sigma","",0.01);
          _n1 = new RooRealVar("n1","",2.);
          _n2 = new RooRealVar("n2","",2.);
          
          _sigmalim = new RooRealConstraint("sigmalim","",*_sigma,0.0002,0.5);
          _meanlim = new RooRealConstraint("meanlim","",*_mean,0.2,2.0);
          _n1lim = new RooRealConstraint("n1lim","",*_n1,1.01,5000.);
          _n2lim = new RooRealConstraint("n2lim","",*_n2,1.01,5000.);
          
          RooConstVar *alpha1 = new RooConstVar("alpha1","",2.0);
          RooConstVar *alpha2 = new RooConstVar("alpha2","",1.0);
          
          _pdf = new RooDoubleCBFast("sigpdf","",*_tgt,*_meanlim,*_sigmalim,*alpha1,*_n1lim,*alpha2,*_n2lim);
          
          //add to RooArgList for proper garbage collection
          _args.addOwned(*_tgt);
          _args.addOwned(*_mean);
          _args.addOwned(*_sigma);
          _args.addOwned(*_n1);
          _args.addOwned(*_n2);
          _args.addOwned(*alpha1);
          _args.addOwned(*alpha2);
          _args.addOwned(*_sigmalim);
          _args.addOwned(*_meanlim);
          _args.addOwned(*_n1lim);
          _args.addOwned(*_n2lim);
          _args.addOwned(*_pdf);        

       }
       else {  
        std::cout << "PhotonAnalysis -- Regression versions 5 and 8 are implemented only!" << std::endl;
        (assert(0)); 
       } 
   }
    // --------------------------------------------------------------------
   if(PADEBUG)
       cout << "InitRealPhotonAnalysis END"<<endl;

   // FIXME book of additional variables

    if(optimizeMVA){
        // Initialize all MVA ---------------------------------------------------//
        l.SetAllMVA();
        // UCSD 
        l.tmvaReaderID_UCSD->BookMVA("Gradient"      ,photonLevelMvaUCSD.c_str()  );
        l.tmvaReader_dipho_UCSD->BookMVA("Gradient"  ,eventLevelMvaUCSD.c_str()   );
        // New ID MVA 
        if( photonLevel2012IDMVA_EB != "" && photonLevel2012IDMVA_EE != "" ) {
            l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevel2012IDMVA_EB.c_str());
            l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevel2012IDMVA_EE.c_str());
        } else {
            assert( run7TeV4Xanalysis );
        }
        // MIT 
        if( photonLevel2011IDMVA_EB != "" && photonLevel2011IDMVA_EE != "" ) {
            l.tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost",photonLevel2011IDMVA_EB.c_str());
            l.tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost",photonLevel2011IDMVA_EE.c_str());
        } else {
            assert( ! run7TeV4Xanalysis );
        }
        l.tmvaReader_dipho_MIT->BookMVA("Gradient"   ,eventLevelMvaMIT.c_str()    );
        // ----------------------------------------------------------------------//    
    }
    
    if (bdtTrainingType == "") {
        bdtTrainingType = bdtTrainingPhilosophy;
    }
    if( run7TeV4Xanalysis ) {
        bdtTrainingType = "Old7TeV";
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::setupEscaleSmearer()
{
    if( splitEscaleSyst ) {
        EnergySmearer::energySmearingParameters::eScaleVector tmp_scale_offset;
        EnergySmearer::energySmearingParameters::phoCatVector tmp_scale_cat;
        readEnergyScaleOffsets(scale_offset_corr_error_file, tmp_scale_offset, tmp_scale_cat,false);
        assert( tmp_scale_offset.size() == 1); assert( ! tmp_scale_cat.empty() );
        
        eScaleCorrPars.categoryType = "Automagic";
        eScaleCorrPars.byRun = false;
        eScaleCorrPars.n_categories = tmp_scale_cat.size();
        eScaleCorrPars.photon_categories = tmp_scale_cat;
        
        eScaleCorrPars.scale_offset = tmp_scale_offset[0].scale_offset;
        eScaleCorrPars.scale_offset_error = tmp_scale_offset[0].scale_offset_error;
        
        eScaleCorrPars.smearing_sigma = tmp_scale_offset[0].scale_offset;
        eScaleCorrPars.smearing_stocastic_sigma = tmp_scale_offset[0].scale_stocastic_offset;
        eScaleCorrPars.smearing_sigma_error = tmp_scale_offset[0].scale_offset_error;
        
        EnergySmearer::energySmearingParameters::phoCatVectorIt icat = tmp_scale_cat.begin();
        for( ; icat != tmp_scale_cat.end(); ++icat ) {
            EnergySmearer * theSmear = new EnergySmearer( eScaleSmearer, EnergySmearer::energySmearingParameters::phoCatVector(1,*icat) );
            theSmear->name( eScaleSmearer->name()+"_"+icat->name );
            theSmear->syst_only(true);
            std::cout << "Uncorrelated single photon category smearer " << theSmear->name() << std::endl;
            photonSmearers_.push_back(theSmear);
            eScaleSmearers_.push_back(theSmear);
        }
        
        eScaleCorrSmearer = new EnergySmearer( eScaleCorrPars );
        eScaleCorrSmearer->name("E_scaleCorr");
        eScaleCorrSmearer->doEnergy(true);
        eScaleCorrSmearer->scaleOrSmear(true);
        eScaleCorrSmearer->syst_only(true);
        photonSmearers_.push_back(eScaleCorrSmearer);
        eScaleSmearers_.push_back(eScaleCorrSmearer);
        
        std::cout << "Uncorrelated single photon category smearer " << eScaleCorrSmearer->name() << std::endl;
        
    } else {
        photonSmearers_.push_back(eScaleSmearer);
        eScaleSmearers_.push_back(eScaleSmearer);
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::setupEscaleSyst(LoopAll &l)
{
    for(std::vector<EnergySmearer*>::iterator ei=eScaleSmearers_.begin(); ei!=eScaleSmearers_.end(); ++ei){
        systPhotonSmearers_.push_back( *ei );
        std::vector<std::string> sys(1,(*ei)->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::addResolSmearer(EnergySmearer * theSmear)
{
    std::cout << "Adding enegy resolution smearer " << theSmear->name() << std::endl;
    EnergySmearerExtrapolation * theExtra = 0;
    if( doStocasticSmearingSyst ) {
        theExtra = new EnergySmearerExtrapolation(theSmear);
        if( ! theExtra->needed() ) { 
            delete theExtra;
            theExtra = 0;
        } else {
            std::cout << "Adding corresponding extrapolation uncertainty term " << theExtra->name() << std::endl;
            photonSmearers_.push_back(theExtra);
        }
    }
    photonSmearers_.push_back(theSmear);
    eResolSmearers_.push_back(std::make_pair(theSmear,theExtra));
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::setupEresolSmearer()
{
    if( splitEresolSyst ) {
        // Use the same format used for the run-dependent energy corrections
        EnergySmearer::energySmearingParameters::eScaleVector tmp_smearing;
        EnergySmearer::energySmearingParameters::phoCatVector tmp_smearing_cat;
        readEnergyScaleOffsets(corr_smearing_file, tmp_smearing, tmp_smearing_cat,false);
        
        // make sure that the scale correction and smearing info is as expected
        assert( tmp_smearing.size() == 1 );
        assert( ! tmp_smearing_cat.empty() );
        
        // copy the read info to the smarer parameters
        eResolCorrPars.categoryType = "Automagic";
        eResolCorrPars.byRun = false;
        eResolCorrPars.n_categories = tmp_smearing_cat.size();
        eResolCorrPars.photon_categories = tmp_smearing_cat;

        eResolCorrPars.scale_offset = tmp_smearing[0].scale_offset;
        eResolCorrPars.scale_offset_error = tmp_smearing[0].scale_offset_error;
        
        eResolCorrPars.smearing_sigma = tmp_smearing[0].scale_offset;
        eResolCorrPars.smearing_sigma_error = tmp_smearing[0].scale_offset_error;
        eResolCorrPars.smearing_stocastic_sigma = tmp_smearing[0].scale_stocastic_offset;
        eResolCorrPars.smearing_stocastic_sigma_error = tmp_smearing[0].scale_stocastic_offset_error;
        eResolCorrPars.smearing_stocastic_pivot = tmp_smearing[0].scale_stocastic_pivot;
        
        EnergySmearer::energySmearingParameters::phoCatVectorIt icat = tmp_smearing_cat.begin();
        for( ; icat != tmp_smearing_cat.end(); ++icat ) {
            EnergySmearer * theSmear = new EnergySmearer( eResolSmearer, EnergySmearer::energySmearingParameters::phoCatVector(1,*icat) );
            theSmear->name( eResolSmearer->name()+"_"+icat->name );
            addResolSmearer(theSmear);
        }
        
        eResolCorrSmearer = new EnergySmearer( eResolCorrPars );
        eResolCorrSmearer->name("E_resCorr");
        eResolCorrSmearer->doEnergy(true);
        eResolCorrSmearer->scaleOrSmear(true);
        eResolCorrSmearer->syst_only(true);
        addResolSmearer(eResolCorrSmearer);
    } else {
        addResolSmearer(eResolSmearer);
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::setupEresolSyst(LoopAll &l)
{
    for(std::vector<std::pair<EnergySmearer*,EnergySmearerExtrapolation*> >::iterator ei=eResolSmearers_.begin(); 
        ei!=eResolSmearers_.end(); ++ei){
        EnergySmearerExtrapolation * extra = ei->second;
        if( extra ) { 
            systPhotonSmearers_.push_back( extra );
            std::vector<std::string> sys(1,extra->name());
            std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
            l.rooContainer->MakeSystematicStudy(sys,sys_t);
        }
        EnergySmearer * sme = ei->first;
        systPhotonSmearers_.push_back( sme );
        std::vector<std::string> sys(1,sme->name());
        std::vector<int> sys_t(1,-1);   // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
}

// ----------------------------------------------------------------------------------------------------
bool PhotonAnalysis::Analysis(LoopAll& l, Int_t jentry)
{
    if(PADEBUG)
        cout << "Analysis START"<<endl;
    pho_presel.clear();

    //remove process ID 18 from gamma+jet to avoid double counting with born+box
    if (l.itype[l.current]==3 && l.process_id==18) return false;

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

    if( pho_acc.size() < 2 || pho_et[ pho_acc[0] ] < presel_scet1 ) return false;

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


        //Fill histograms according to diphoton or single photon category, as appropriate

        int dipho_category = l.DiphotonCategory(diphoton_index.first, diphoton_index.second, Higgs.Pt(), Higgs.Pt()/Higgs.M(), 2, 2, 2);
        int leadpho_category = l.PhotonCategory(diphoton_index.first, 2, 2);
        int subleadpho_category = l.PhotonCategory(diphoton_index.second, 2, 2);


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

    return true;
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
                                   *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])),
                                   ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(),
                                   energyCorrected[ipho],
                                   l.pho_isEB[ipho],
                                   l.pho_r9[ipho],
                                   ipho,
                                   false,
                                   (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
                                   );
        float pweight = 1.;
        float sweight = 1.;
        float eta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[ipho]))->Eta());

        if( doEcorrectionSmear )  {
            eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
        }
        /// FIXME deterministic smearing on MC photons
        if( cur_type == 0 ) {          // correct energy scale in data
            float ebefore = phoInfo.energy();
            eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
            pweight *= sweight;
        }
        // apply mc-derived photon corrections, to data and MC alike
        corrected_pho_energy[ipho] = phoInfo.energy();
    }

    for(int ipho=0; ipho<l.pho_n; ++ipho) {

        l.pho_genmatched[ipho]=GenMatchedPhoton( l, ipho);

        // match all photons in the original tree with the conversions from the merged collection and save the indices
        int iConv  =l.matchPhotonToConversion(ipho);
        if ( iConv>=0 )
            (*l.pho_matchingConv).push_back(l.matchPhotonToConversion(ipho));
        else
            (*l.pho_matchingConv).push_back(-1);

        TLorentzVector p4 = l.get_pho_p4(ipho,0,&corrected_pho_energy[0]);
        float eta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[ipho]))->Eta());
        // photon et wrt 0,0,0
        float et = p4.Pt();
        pho_et.push_back(et);

	if( (eta>1.4442 && eta<1.566) || eta > 2.5 ) {
            continue;
        }
        pho_acc.push_back(ipho);

        pho_presel.push_back(ipho);
    }

    std::sort(pho_acc.begin(),pho_acc.end(),
              SimpleSorter<float,std::greater<float> >(&pho_et[0]));
    std::sort(pho_presel.begin(),pho_presel.end(),
              SimpleSorter<float,std::greater<float> >(&pho_et[0]));

    if( pho_presel.size() > 1 ) {
        for(size_t ipho=0; ipho<pho_presel.size()-1; ++ipho ) {
            assert( pho_et[pho_presel[ipho]] >= pho_et[pho_presel[ipho+1]] );
        }
    }
    if( pho_acc.size()>1 ) {
        for(size_t ipho=0; ipho<pho_acc.size()-1; ++ipho ) {
            assert( pho_et[pho_acc[ipho]] >= pho_et[pho_acc[ipho+1]] );
        }
    }

}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
    if(PADEBUG)
	cout<<"myFillReduceVar START"<<endl;

    // Run on-the-fly regression at Reduction Step
    if( l.typerun == LoopAll::kReduce ) {
        GetRegressionCorrections(l);  // need to pass LoopAll
    }
    PreselectPhotons(l,jentry);

    if(PADEBUG)
        cout<<"myFillReduceVar END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::postProcessJets(LoopAll & l, int vtx)
{
    int minv = 0, maxv = l.vtx_std_n;
    if( l.typerun == l.kFill && l.version > 14 && maxv >= l.jet_algoPF1_nvtx ) {
	maxv = l.jet_algoPF1_nvtx-1;
    }
    if( vtx!=-1 ){
	minv = vtx;
	maxv = vtx+1;
    }
    if( vtx == -1 || vtx == 0 ) {
	for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) {
	    if( recorrectJets ) {
		jetHandler_->recomputeJec(ijet, true);
	    }
	    if( emulateJetResponse && l.itype[l.current] != 0 ) {
		jetHandler_->emulateJetResponse(ijet);
	    }
	    if( applyJer ) {
		jetHandler_->applyJerUncertainty(ijet, jerShift);
	    }
	    if( applyJecUnc ) {
		jetHandler_->applyJecUncertainty(ijet, jecShift);
	    }
	}
    }
    for(int ivtx=minv;ivtx<maxv; ++ivtx) {
	for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) {
	    if( recomputeBetas || (l.typerun != l.kFill && l.version > 14 && ivtx >= l.jet_algoPF1_nvtx) ) {
		/// std::cout << "recomputeBetas " << ivtx << " " << l.jet_algoPF1_nvtx << std::endl;
		jetHandler_->computeBetas(ijet, ivtx);
	    }
	    if( rerunJetMva ) {
		jetHandler_->computeMva(ijet, ivtx);
	    } else if ( recomputeJetWp ) {
		jetHandler_->computeWp(ijet, ivtx);
	    }
	}
	if( ivtx >= l.jet_algoPF1_nvtx && l.typerun != l.kFill ) {
	    l.jet_algoPF1_nvtx = ivtx+1;
	}
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::switchJetIdVertex(LoopAll &l, int ivtx)
{
    /// if( l.version > 14 ) { l.jet_algoPF1_nvtx = 10; }
    if( l.jet_algoPF1_n > 0 && l.jet_algoPF1_nvtx < (*l.jet_algoPF1_betaStarClassic_ext)[0].size() ) {
	l.jet_algoPF1_nvtx = (*l.jet_algoPF1_betaStarClassic_ext)[0].size();
    }
    if( l.typerun == l.kFill && l.version > 14 && ivtx >= l.jet_algoPF1_nvtx ) {
	std::cout << "WARNING choosen vertex beyond " << l.jet_algoPF1_nvtx << " and jet ID was not computed. Falling back to vertex 0." << std::endl;
	ivtx = 0;
    }
    //cout << " ivtx " << ivtx << " jalgonvtx " << l.jet_algoPF1_nvtx << endl;

    for(int ii=0; ii<l.jet_algoPF1_n; ++ii) {
	l.jet_algoPF1_beta[ii]              = (*l.jet_algoPF1_beta_ext)[ii][ivtx];
        l.jet_algoPF1_betaStar[ii]          = (*l.jet_algoPF1_betaStar_ext)[ii][ivtx];
        l.jet_algoPF1_betaStarClassic[ii]   = (*l.jet_algoPF1_betaStarClassic_ext)[ii][ivtx];
        l.jet_algoPF1_simple_mva[ii]        = (*l.jet_algoPF1_simple_mva_ext)[ii][ivtx];
        l.jet_algoPF1_full_mva[ii]          = (*l.jet_algoPF1_full_mva_ext)[ii][ivtx];
        l.jet_algoPF1_simple_wp_level[ii]   = (*l.jet_algoPF1_simple_wp_level_ext)[ii][ivtx];
        l.jet_algoPF1_full_wp_level[ii]     = (*l.jet_algoPF1_full_wp_level_ext)[ii][ivtx];
        l.jet_algoPF1_cutbased_wp_level[ii] = (*l.jet_algoPF1_cutbased_wp_level_ext)[ii][ivtx];
    }
}


// ----------------------------------------------------------------------------------------------------
bool PhotonAnalysis::SelectEventsReduction(LoopAll& l, int jentry)
{

    if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
    // require at least two reconstructed photons to store the event

    if( pho_acc.size() < 2 ) { return false; }

    vtxAna_.clear();
    l.vtx_std_ranked_list->clear();
    l.dipho_vtx_std_sel->clear();
    l.vtx_std_ranked_list->clear();
    l.vtx_std_evt_mva->clear();
    l.vtx_std_sel=0;
    float maxSumPt = 0.;
    l.dipho_n = 0;
    bool oneKinSelected = false;

    // fill ID variables
    if( forcedRho >= 0. ) {
        l.rho = forcedRho;
    } else {
        l.rho = l.rho_algo1;
    }
    
    l.FillCICInputs();
    if(reComputeCiCPF) { l.FillCICPFInputs(); }
    l.FillCIC();
    l.FillMuonGsfTracks();

    if(l.itype[l.current]<0) {
        bool foundHiggs=FindHiggsObjects(l);
        if(PADEBUG)  cout << " foundHiggs? "<<foundHiggs<<std::endl;
    } else {
        SetNullHiggs(l);
    }
    /// Jet matching
    // pfJets ak5
    l.doJetMatching(*l.jet_algoPF1_p4,*l.genjet_algo1_p4,l.jet_algoPF1_genMatched,l.jet_algoPF1_vbfMatched,l.jet_algoPF1_bgenMatched,l.jet_algoPF1_cgenMatched,l.jet_algoPF1_lgenMatched,l.jet_algoPF1_genPt,l.jet_algoPF1_genDr);
    // pfJets ak7
    //l.doJetMatching(*l.jet_algoPF2_p4,*l.genjet_algo2_p4,l.jet_algoPF2_genMatched,l.jet_algoPF2_vbfMatched,l.jet_algoPF2_genPt,l.jet_algoPF2_genDr);
    // CHS ak5
    l.doJetMatching(*l.jet_algoPF3_p4,*l.genjet_algo1_p4,l.jet_algoPF3_genMatched,l.jet_algoPF3_vbfMatched,l.jet_algoPF3_bgenMatched,l.jet_algoPF3_cgenMatched,l.jet_algoPF3_lgenMatched,l.jet_algoPF3_genPt,l.jet_algoPF3_genDr);    

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

void PhotonAnalysis::MetCorrections2012(LoopAll& l)
{
    //shift met (reduction step)
    //in mc smearing should be applied first and then shifting
    //both these are performed at analysis step
    TLorentzVector unpfMET;
    unpfMET.SetPxPyPzE (l.met_pfmet*cos(l.met_phi_pfmet),
           l.met_pfmet*sin(l.met_phi_pfmet),
           0,
           sqrt(l.met_pfmet*cos(l.met_phi_pfmet) * l.met_pfmet*cos(l.met_phi_pfmet)
           + l.met_pfmet*sin(l.met_phi_pfmet) * l.met_pfmet*sin(l.met_phi_pfmet)));

    bool isMC = l.itype[l.current]!=0;

    TLorentzVector shiftMET_corr = l.shiftMet(&unpfMET,isMC);
    l.shiftMET_pt = shiftMET_corr.Pt();
    l.shiftMET_phi = shiftMET_corr.Phi();
    l.shiftMET_eta = shiftMET_corr.Eta();
    l.shiftMET_e = shiftMET_corr.Energy();
}


void PhotonAnalysis::MetCorrections2012_Simple(LoopAll& l,TLorentzVector lead_p4 ,TLorentzVector sublead_p4)
{
    // mc: scaling and shifting, data: scaling (analysis step)
    TLorentzVector unpfMET;
    unpfMET.SetPxPyPzE (l.met_pfmet*cos(l.met_phi_pfmet),
           l.met_pfmet*sin(l.met_phi_pfmet),
           0,
           sqrt(l.met_pfmet*cos(l.met_phi_pfmet) * l.met_pfmet*cos(l.met_phi_pfmet)
           + l.met_pfmet*sin(l.met_phi_pfmet) * l.met_pfmet*sin(l.met_phi_pfmet)));

     bool isMC = l.itype[l.current]!=0;

     //take shifted met for data
     TLorentzVector shiftedMET;
     double shiftedMETpt = l.shiftMET_pt;
     double shiftedMETe = l.shiftMET_e;
     double shiftedMETeta = l.shiftMET_eta;
     double shiftedMETphi = l.shiftMET_phi;

     shiftedMET.SetPtEtaPhiE(shiftedMETpt,shiftedMETeta,shiftedMETphi,shiftedMETe);

     if (isMC) {
       //smear raw met for mc
       TLorentzVector smearMET_corr = l.correctMet_Simple( lead_p4, sublead_p4 , &unpfMET, true, false);
       l.smearMET_pt = smearMET_corr.Pt();
       l.smearMET_phi = smearMET_corr.Phi();
       //shift smeared met for mc
       TLorentzVector shiftsmearMET_corr = l.shiftMet(&smearMET_corr,isMC);
       l.shiftsmearMET_pt = shiftsmearMET_corr.Pt();
       l.shiftsmearMET_phi = shiftsmearMET_corr.Phi();
       l.correctedpfMET = l.shiftsmearMET_pt;
       l.correctedpfMET_phi = l.shiftsmearMET_phi;
     } else {
       //scale shifted met for data
       TLorentzVector shiftscaleMET_corr = l.correctMet_Simple( lead_p4, sublead_p4 , &shiftedMET, false , true);
       l.shiftscaleMET_pt = shiftscaleMET_corr.Pt();
       l.shiftscaleMET_phi = shiftscaleMET_corr.Phi();
       l.correctedpfMET = l.shiftscaleMET_pt;
       l.correctedpfMET_phi = l.shiftscaleMET_phi;
     }
}


bool PhotonAnalysis::SkimEvents(LoopAll& l, int jentry)
{
    static TH1F * promptFakeFractions = 0;
    static TH1F * promptMotherStatus = 0;
    static TH1F * fakeMotherStatus = 0;
    if( run7TeV4Xanalysis ) { l.version=12; }
    
    
    l.b_pho_n->GetEntry(jentry);
    if( l.pho_n < 2 ) {
        return false;
    }

    if( skimOnDiphoN && l.typerun == l.kFill ) {
	l.b_dipho_n->GetEntry(jentry);
	if( l.dipho_n < 1 ) {
	    return false;
	}
    }

    // do not run trigger selection on MC
    int filetype = l.itype[l.current];
    bool skipTrigger = !doTriggerSelection || ( filetype != 0 && !l.runZeeValidation ) || triggerSelections.empty() || (l.sqrtS == 7 && filetype != 0);
    
    if( ! skipTrigger ) {
        // get the trigger selection for this run
        l.b_run->GetEntry(jentry);
        std::vector<TriggerSelection>::iterator isel = find(triggerSelections.begin(), triggerSelections.end(), l.run );
        if(isel == triggerSelections.end() ) {
            std::cerr << "No trigger selection for run " << l.run << "defined" << std::endl;
            return true;
        }

	// get the trigger data
        if( l.version < 13 ) {
            l.b_hlt1_bit->GetEntry(jentry);
            l.b_hlt_path_names_HLT1->GetEntry(jentry);
            if( !  isel->pass( *(l.hlt_path_names_HLT1), *(l.hlt1_bit) ) ) {
                return false;
            }
        } else {
            l.b_hlt_bit->GetEntry(jentry);
            l.b_hlt_path_names_HLT->GetEntry(jentry);
            if( !  isel->pass( *(l.hlt_path_names_HLT), *(l.hlt_bit) ) ) {
                return false;
            }
        }
        //l.countersred[trigCounter_]++;
    }

    if( l.typerun == l.kReduce || l.typerun == l.kFillReduce ) {
        //// if( filetype == 2 ) { // photon+jet
        ////    l.b_process_id->GetEntry(jentry);
        ////    if( l.process_id == 18 ) {
        ////        return false;
        ////    }
        //// }

        if(selectprocess){
            if(processtoselect!=l.process_id){
                return false;
            }
        }

        if( filetype != 0 && ! (keepPP && keepPF && keepFF) ) {
            if( promptFakeFractions == 0 ) {
                promptFakeFractions = new TH1F("promptFakeFractions","promptFakeFractions",3,-0.5,2.5);
                promptMotherStatus = new TH1F("promptMotherStatus","promptMotherStatus",20,-0.5,20);
                fakeMotherStatus = new TH1F("fakeMotherStatus","fakeMotherStatus",20,-0.5,20);
                l.AddGlobalHisto(promptFakeFractions);
                l.AddGlobalHisto(promptMotherStatus);
                l.AddGlobalHisto(fakeMotherStatus);
            }
            
            l.b_weight->GetEntry(jentry);
            l.b_gp_n->GetEntry(jentry);
            l.b_gp_mother->GetEntry(jentry);
            l.b_gp_status->GetEntry(jentry);
            l.b_gp_pdgid->GetEntry(jentry);
            l.b_gp_p4->GetEntry(jentry);

            std::vector<TLorentzVector *> gen_pho;
            
            int np = 0;
            for(int ip=0;ip<l.gp_n;++ip) {
                if( l.gp_status[ip] != 1 || l.gp_pdgid[ip] != 22 ) {
                    continue;
                }
                TLorentzVector * p4 = (TLorentzVector*) l.gp_p4->At(ip);
                if( l.gp_mother[ip] < 0 || p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
                bool duplicate = false;
                for(size_t ii=0; ii<gen_pho.size(); ++ii) {
                    if( p4->DeltaR(*p4) < 0.05 ) { 
                        duplicate = true; 
                        break;
                    }
                }
                if( duplicate ) { continue; }
                int mother_id = abs( l.gp_pdgid[ l.gp_mother[ip] ] );
                if( mother_id <= 25 ) { 
                    ++np; 
                    promptMotherStatus->Fill((float)l.gp_status[l.gp_mother[ip]],l.weight);
                } else {
                    fakeMotherStatus->Fill((float)l.gp_status[l.gp_mother[ip]],l.weight);
                }
                if( np >= 2 ) { break; }
            }
            /// std::cout << "N prompt photons: " << np << std::endl;
            promptFakeFractions->Fill((float)np,l.weight);
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
    l.dipho_vtx_std_sel =  new std::vector<int>();

    if( outputTree ) {

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
	l.Branch_jet_algoPF1_bgenMatched(outputTree);
	l.Branch_jet_algoPF1_cgenMatched(outputTree);
	l.Branch_jet_algoPF1_lgenMatched(outputTree);
	l.Branch_jet_algoPF1_vbfMatched(outputTree);
	l.Branch_jet_algoPF1_genPt(outputTree);
	l.Branch_jet_algoPF1_genDr(outputTree);

	//l.Branch_jet_algoPF2_genMatched(outputTree);
	//l.Branch_jet_algoPF2_vbfMatched(outputTree);
	//l.Branch_jet_algoPF2_genPt(outputTree);
	//l.Branch_jet_algoPF2_genDr(outputTree);

	l.Branch_jet_algoPF3_genMatched(outputTree);
	l.Branch_jet_algoPF3_bgenMatched(outputTree);
	l.Branch_jet_algoPF3_cgenMatched(outputTree);
	l.Branch_jet_algoPF3_lgenMatched(outputTree);
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

void PhotonAnalysis::ResetAnalysis(){
}

void PhotonAnalysis::SetNullHiggs(LoopAll& l){

    ((TLorentzVector *)l.gh_higgs_p4->At(0))->SetXYZT(0,0,0,0);

    ((TLorentzVector *)l.gh_pho1_p4->At(0))->SetXYZT(0,0,0,0);
    ((TLorentzVector *)l.gh_pho2_p4->At(0))->SetXYZT(0,0,0,0);

    ((TLorentzVector *)l.gh_vbfq1_p4->At(0))->SetXYZT(0,0,0,0);
    ((TLorentzVector *)l.gh_vbfq2_p4->At(0))->SetXYZT(0,0,0,0);
    l.gh_vbfq1_pdgid=0;
    l.gh_vbfq2_pdgid=0;

    ((TLorentzVector *)l.gh_vh1_p4->At(0))->SetXYZT(0,0,0,0);
    ((TLorentzVector *)l.gh_vh2_p4->At(0))->SetXYZT(0,0,0,0);
    l.gh_vh_pdgid=0;
    l.gh_vh1_pdgid=0;
    l.gh_vh2_pdgid=0;

}


bool PhotonAnalysis::FindHiggsObjects(LoopAll& l){

    int higgsind=-1;
    int mc1=-1;
    int mc2=-1;
    int i1=-1;
    int i2=-1;

    l.FindMCHiggsPhotons( higgsind,  mc1,  mc2,  i1,  i2 );


    if(higgsind!=-1) {
        TLorentzVector * TheHiggs = (TLorentzVector *) l.gp_p4->At(higgsind);
        ((TLorentzVector *)l.gh_higgs_p4->At(0))->SetXYZT(TheHiggs->Px(),TheHiggs->Py(),TheHiggs->Pz(),TheHiggs->E());
    } else { ((TLorentzVector *)l.gh_higgs_p4->At(0))->SetXYZT(0,0,0,0); }

    l.gh_gen2reco1=i1;
    l.gh_gen2reco2=i2;


    if(mc1!=-1) {
        TLorentzVector * mcpho1 = (TLorentzVector *) l.gp_p4->At(mc1);
        ((TLorentzVector *)l.gh_pho1_p4->At(0))->SetXYZT(mcpho1->Px(),mcpho1->Py(),mcpho1->Pz(),mcpho1->E());
    } else { ((TLorentzVector *)l.gh_pho1_p4->At(0))->SetXYZT(0,0,0,0); }

    if(mc2!=-1) {
        TLorentzVector * mcpho2 = (TLorentzVector *) l.gp_p4->At(mc2);
        ((TLorentzVector *)l.gh_pho2_p4->At(0))->SetXYZT(mcpho2->Px(),mcpho2->Py(),mcpho2->Pz(),mcpho2->E());
    } else { ((TLorentzVector *)l.gh_pho2_p4->At(0))->SetXYZT(0,0,0,0); }


    int vbfq1=-100;
    int vbfq2=-100;

    int vh=-100;
    int vh1=-100;
    int vh2=-100;

    if(higgsind!=-1){
        l.FindMCVBF(higgsind,vbfq1,vbfq2);
        l.FindMCVH(higgsind,vh,vh1,vh2);

        //std::cout<<"higgsind vbfq1 vbfq2 vh gp_pdgid[vh] vh1 vh2  "<<higgsind<<"  "<<vbfq1<<"  "<<vbfq2<<"  "<<vh<<"  "<<l.gp_pdgid[vh]<<"  "<<vh1<<"  "<<vh2<<std::endl;
    }


    if(vh==-100) l.gh_vh_pdgid=-10000;
    else l.gh_vh_pdgid=l.gp_pdgid[vh];

    if(vh1==-100) l.gh_vh1_pdgid=-10000;
    else l.gh_vh1_pdgid=l.gp_pdgid[vh1];

    if(vh2==-100) l.gh_vh2_pdgid=-10000;
    else l.gh_vh2_pdgid=l.gp_pdgid[vh2];

    if(vh1!=-100) {
        TLorentzVector * mcvh1 = (TLorentzVector *) l.gp_p4->At(vh1);
        ((TLorentzVector *)l.gh_vh1_p4->At(0))->SetXYZT(mcvh1->Px(),mcvh1->Py(),mcvh1->Pz(),mcvh1->E());
    } else { ((TLorentzVector *)l.gh_vh1_p4->At(0))->SetXYZT(0,0,0,0); }

    if(vh2!=-100) {
        TLorentzVector * mcvh2 = (TLorentzVector *) l.gp_p4->At(vh2);
        ((TLorentzVector *)l.gh_vh2_p4->At(0))->SetXYZT(mcvh2->Px(),mcvh2->Py(),mcvh2->Pz(),mcvh2->E());
    } else { ((TLorentzVector *)l.gh_vh2_p4->At(0))->SetXYZT(0,0,0,0); }



    if(vbfq1==-100) l.gh_vbfq1_pdgid=0;
    else l.gh_vbfq1_pdgid=l.gp_pdgid[vbfq1];

    if(vbfq2==-100) l.gh_vbfq2_pdgid=0;
    else l.gh_vbfq2_pdgid=l.gp_pdgid[vbfq2];

    if(vbfq1!=-100) {
        TLorentzVector * mcvbfq1 = (TLorentzVector *) l.gp_p4->At(vbfq1);
        ((TLorentzVector *)l.gh_vbfq1_p4->At(0))->SetXYZT(mcvbfq1->Px(),mcvbfq1->Py(),mcvbfq1->Pz(),mcvbfq1->E());
    } else { ((TLorentzVector *)l.gh_vbfq1_p4->At(0))->SetXYZT(0,0,0,0); }

    if(vbfq2!=-100) {
        TLorentzVector * mcvbfq2 = (TLorentzVector *) l.gp_p4->At(vbfq2);
        ((TLorentzVector *)l.gh_vbfq2_p4->At(0))->SetXYZT(mcvbfq2->Px(),mcvbfq2->Py(),mcvbfq2->Pz(),mcvbfq2->E());
    } else { ((TLorentzVector *)l.gh_vbfq2_p4->At(0))->SetXYZT(0,0,0,0); }

    return (higgsind != -1);

}



Bool_t PhotonAnalysis::GenMatchedPhoton(LoopAll& l, int ipho){
    Bool_t is_prompt = false;
    TLorentzVector* phop4 = (TLorentzVector*) l.pho_p4->At(ipho);

    for(int ip=0;ip<l.gp_n;++ip) {
        if( l.gp_status[ip] != 1 || l.gp_pdgid[ip] != 22 ) {
            continue;
        }
        TLorentzVector * p4 = (TLorentzVector*) l.gp_p4->At(ip);
        if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
        int mother_id = abs( l.gp_pdgid[ l.gp_mother[ip] ] );
        if( mother_id <= 25 ) {
            float dr = phop4->DeltaR(*p4);
            if (dr<0.3 && fabs((p4->Pt()-phop4->Pt())/p4->Pt()) < 0.5) {
                is_prompt = true;
                break;
            }
        }
    }
    return is_prompt;
}


bool PhotonAnalysis::ClassicCatsNm1Plots(LoopAll& l, int diphoton_nm1_id, float* smeared_pho_energy, float eventweight, float myweight){

    bool pass = false;
    if(diphoton_nm1_id==-1) return pass;

    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_nm1_id], l.dipho_vtxind[diphoton_nm1_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_nm1_id], l.dipho_vtxind[diphoton_nm1_id], &smeared_pho_energy[0]);
    TLorentzVector diphoton = lead_p4+sublead_p4;

    int photon_category   = l.PhotonCategory(l.dipho_subleadind[diphoton_nm1_id],2,2);
    sublead_r9            = l.pho_r9[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_trkisooet     = (*l.pho_tkiso_recvtx_030_002_0000_10_01)[l.dipho_subleadind[diphoton_nm1_id]][l.dipho_vtxind[diphoton_nm1_id]]*50/sublead_p4.Et();
    sublead_isoOverEt     = (l.pho_hcalsumetconedr04[l.dipho_subleadind[diphoton_nm1_id]]
                             +  l.pho_ecalsumetconedr03[l.dipho_subleadind[diphoton_nm1_id]]
                             +  (*l.pho_tkiso_recvtx_030_002_0000_10_01)[l.dipho_subleadind[diphoton_nm1_id]][l.dipho_vtxind[diphoton_nm1_id]]
                             - 0.17*l.rho)*50/sublead_p4.Et();
    sublead_badisoOverEt  = (l.pho_hcalsumetconedr04[l.dipho_subleadind[diphoton_nm1_id]]
                             +  l.pho_ecalsumetconedr04[l.dipho_subleadind[diphoton_nm1_id]]
                             +  l.pho_tkiso_badvtx_040_002_0000_10_01[l.dipho_subleadind[diphoton_nm1_id]]
                             - 0.52*l.rho)*50/sublead_p4.Et();

    sublead_sieie         = l.pho_sieie[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_drtotk        = l.pho_drtotk_25_99[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_hovere        = l.pho_hoe[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_mgg           = diphoton.M();

    int applyCutsType = 15 + photon_category;

    pass = l.ApplyCutsFill(0,applyCutsType, eventweight, myweight);

    return pass;
}


bool PhotonAnalysis::ElectronTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphotonVHlep_id==-1) return tag;
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);

    int elInd = l.ElectronSelection(lead_p4, sublead_p4, l.dipho_vtxind[diphotonVHlep_id]);
    if(elInd!=-1) tag = true;

    return tag;
}

bool PhotonAnalysis::ElectronTag2012(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, ofstream& lep_sync, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphotonVHlep_id==-1) return tag;
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);

    int elInd = l.ElectronSelection2012(lead_p4, sublead_p4, l.dipho_vtxind[diphotonVHlep_id]);
    if(elInd!=-1) tag = true;

        TLorentzVector dipho_p4 = lead_p4+sublead_p4;
        lep_sync<<"run="<<l.run<<"\t";
        lep_sync<<"lumis"<<l.lumis<<"\t";
        lep_sync<<"event="<<(unsigned int) l.event<<"\t";
        lep_sync<<"pt1="<<lead_p4.Pt()<<"\t";
        lep_sync<<"eta1="<<lead_p4.Eta()<<"\t";
        lep_sync<<"pt2="<<sublead_p4.Pt()<<"\t";
        lep_sync<<"eta2="<<sublead_p4.Eta()<<"\t";
        lep_sync<<"ptgg="<<dipho_p4.Pt()<<"\t";
    if(tag){
        TLorentzVector* el_p4 = (TLorentzVector*) l.el_std_p4->At(elInd);
        lep_sync<<"elpt="<<el_p4->Pt()<<"\t";
        lep_sync<<"eleta="<<el_p4->Eta()<<"  ELTAG\n";
    } else {
        lep_sync<<"ELNOTAG\n";
    }

    return tag;
}



bool PhotonAnalysis::ElectronTag2012B(LoopAll& l, int& diphotonVHlep_id, int& el_ind, int& elVtx, int& el_cat, float* smeared_pho_energy, ofstream& lep_sync, bool mvaselection, float phoidMvaCut, float eventweight, std::vector<float> smeared_pho_weight, bool fillHist, bool vetodipho, bool kinonly){
    bool tag = false;
    float elptcut=20;
    bool localdebug=false;

    el_ind=l.ElectronSelectionMVA2012(elptcut);

    if(el_ind!=-1) {
        if(localdebug) cout<<"in ElectronTag2012B and selected "<<el_ind<<endl;
        TLorentzVector* myel = (TLorentzVector*) l.el_std_p4->At(el_ind);
        TLorentzVector* myelsc = (TLorentzVector*) l.el_std_sc->At(el_ind);
        elVtx=l.FindElectronVertex(el_ind);
        if(localdebug) cout<<"in ElectronTag2012B and selected vtx "<<elVtx<<endl;

        float drtoveto = 0.5;
        std::vector<bool> veto_indices;
        veto_indices.clear();
        l.PhotonsToVeto(myelsc, drtoveto, veto_indices, true);
        for(int iveto=0; iveto<veto_indices.size(); iveto++){
            if(localdebug) cout<<"veto ipho "<<veto_indices[iveto]<<" "<<iveto<<endl;
        }

        if(mvaselection) {
            diphotonVHlep_id = l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVHlepCut,subleadEtVHlepCut,phoidMvaCut,
                applyPtoverM, &smeared_pho_energy[0], vetodipho, kinonly, diphobdt_output_Cut_VHLep, -1, false, veto_indices);
            if(localdebug) cout<<"diphotonVHlep_id "<<diphotonVHlep_id<<endl;
        } else {
            diphotonVHlep_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHlepCut,subleadEtVHlepCut, 4,
                applyPtoverM, &smeared_pho_energy[0], true, -1, veto_indices);
        }

        if(diphotonVHlep_id!=-1){
            TLorentzVector lead_p4      = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id],    l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
            TLorentzVector sublead_p4   = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
            TLorentzVector dipho_p4 = lead_p4 + sublead_p4;
            float mass = dipho_p4.M();

            // need to check again for d0 and dZ (couldn't before because we didn't have the vertex)
            if(l.ElectronMVACuts(el_ind, elVtx)){
                std::string label("noleppho_nomva");
                if(mass>=100 && mass<180 && fillHist){
                    if(smeared_pho_weight.size()!=0){
                        if(localdebug) cout<<"ElectronTag2012B eventweight passed in "<<eventweight<<std::endl;
                        if(localdebug) cout<<"meared_pho_weight l.dipho_leadind[diphotonVHlep_id] diphotonVHlep_id "<<smeared_pho_weight[l.dipho_leadind[diphotonVHlep_id]]<<" "<<smeared_pho_weight[l.dipho_subleadind[diphotonVHlep_id]]<<" "<<l.dipho_leadind[diphotonVHlep_id]<<" "<<l.dipho_subleadind[diphotonVHlep_id]<<" "<<diphotonVHlep_id<<std::endl;
                        eventweight*=(smeared_pho_weight[l.dipho_leadind[diphotonVHlep_id]] * smeared_pho_weight[l.dipho_subleadind[diphotonVHlep_id]]);
                        if(localdebug) cout<<"ElectronTag2012B eventweight*smearedphoweight "<<eventweight<<std::endl;
                    }
                    int cur_type = l.itype[l.current];
                    ControlPlotsElectronTag2012B(l, lead_p4, sublead_p4, el_ind, 0., eventweight, label);
                }

                if(l.ElectronPhotonCuts2012B(lead_p4, sublead_p4, *myel, includeVHlepPlusMet)){ 
                    tag=true;
                    el_cat=(int)(abs(lead_p4.Eta())>1.5 || abs(sublead_p4.Eta())>1.5);
                    if(localdebug) cout<<"pass ElectronPhotonCuts2012B, el_cat "<<el_cat<<endl;
                } else {
                    diphotonVHlep_id=-1;
                    if(localdebug) cout<<"fail ElectronPhotonCuts2012B"<<endl;
                }
            }

            if(tag){
                lep_sync<<"run="<<l.run<<"\t";
                lep_sync<<"lumis"<<l.lumis<<"\t";
                lep_sync<<"event="<<(unsigned int) l.event<<"\t";
                lep_sync<<"pt1="<<lead_p4.Pt()<<"\t";
                lep_sync<<"eta1="<<lead_p4.Eta()<<"\t";
                lep_sync<<"pt2="<<sublead_p4.Pt()<<"\t";
                lep_sync<<"eta2="<<sublead_p4.Eta()<<"\t";
                lep_sync<<"ptgg="<<dipho_p4.Pt()<<"\t";
                lep_sync<<"elpt="<<myel->Pt()<<"\t";
                lep_sync<<"eleta="<<myel->Eta()<<"  ELTAG\n";
            }
        }
    }

    return tag;
}

void PhotonAnalysis::ControlPlotsElectronTag2012B(LoopAll& l, TLorentzVector lead_p4, TLorentzVector sublead_p4, int el_ind, float bdtoutput, float evweight, std::string label){
    if(el_ind<0) {
        std::cout<<"el_ind is "<<el_ind<<std::endl;
        std::cout<<"leaving ControlPlotsElectronTag2012B"<<std::endl;
        return;
    }
    TLorentzVector* myel = (TLorentzVector*) l.el_std_p4->At(el_ind);
    TLorentzVector* myelsc = (TLorentzVector*) l.el_std_sc->At(el_ind);
    int elVtx=l.FindElectronVertex(el_ind);

    double Aeff=0.;
    float thiseta = fabs(myelsc->Eta());
    if(thiseta<1.0)                   Aeff=0.135;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.168;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.068;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.116;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.162;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.241;
    if(thiseta>=2.4)                  Aeff=0.23;
    float thisiso=l.el_std_pfiso_charged[el_ind]+std::max(l.el_std_pfiso_neutral[el_ind]+l.el_std_pfiso_photon[el_ind]-l.rho*Aeff,0.);
    TLorentzVector elead_p4 = *myel + lead_p4;
    TLorentzVector esub_p4 = *myel + sublead_p4;

    int el_cat=(int)(abs(lead_p4.Eta())>1.5 || abs(sublead_p4.Eta())>1.5);
    l.FillHist(Form("ElectronTag_elpt_%s",label.c_str()),          el_cat, myel->Pt(), evweight);
    l.FillHist(Form("ElectronTag_eleta_%s",label.c_str()),         el_cat, myelsc->Eta(), evweight);
    l.FillHist(Form("ElectronTag_elPhi_%s",label.c_str()),         el_cat, myelsc->Phi(), evweight);
    l.FillHist(Form("ElectronTag_elmva_%s",label.c_str()),         el_cat, l.el_std_mva_nontrig[el_ind], evweight);
    l.FillHist(Form("ElectronTag_elisoopt_%s",label.c_str()),      el_cat, thisiso/myel->Pt(), evweight);
    l.FillHist(Form("ElectronTag_missinghits_%s",label.c_str()),   el_cat, l.el_std_hp_expin[el_ind], evweight);
    l.FillHist(Form("ElectronTag_passconvveto_%s",label.c_str()),  el_cat, l.el_std_conv[el_ind], evweight);
    l.FillHist(Form("ElectronTag_drellead_%s",label.c_str()),      el_cat, myel->DeltaR(lead_p4), evweight);
    l.FillHist(Form("ElectronTag_drelsub_%s",label.c_str()),       el_cat, myel->DeltaR(sublead_p4), evweight);
    l.FillHist(Form("ElectronTag_d0_%s",label.c_str()),            el_cat, l.el_std_D0Vtx[el_ind][elVtx], evweight);
    l.FillHist(Form("ElectronTag_dZ_%s",label.c_str()),            el_cat, l.el_std_DZVtx[el_ind][elVtx], evweight);
    l.FillHist(Form("ElectronTag_Melsub_%s",label.c_str()),        el_cat, esub_p4.M(), evweight);
    l.FillHist(Form("ElectronTag_Mellead_%s",label.c_str()),       el_cat, elead_p4.M(), evweight);
    l.FillHist(Form("ElectronTag_dMelsub_%s",label.c_str()),       el_cat, abs(esub_p4.M()-91.2), evweight);
    l.FillHist(Form("ElectronTag_dMellead_%s",label.c_str()),      el_cat, abs(elead_p4.M()-91.2), evweight);
    l.FillHist(Form("ElectronTag_mindMel_%s",label.c_str()),       el_cat, (float)min(abs(esub_p4.M()-91.2),abs(elead_p4.M()-91.2)), evweight);
    l.FillHist(Form("ElectronTag_diphomva_%s",label.c_str()),      el_cat, bdtoutput, evweight);

}

void PhotonAnalysis::ControlPlotsMetTag2012B(LoopAll& l, TLorentzVector lead_p4, TLorentzVector sublead_p4, float bdtoutput, float evweight, std::string label){
    TLorentzVector myMet = l.METCorrection2012B(lead_p4, sublead_p4, moriond2013MetCorrection);
    float corrMet    = myMet.Pt();
    float corrMetPhi = myMet.Phi();

    //    std::cout << "***corrMet " << corrMet << " corrPhi " << corrMetPhi << std::endl;

    int met_cat=(int)(abs(lead_p4.Eta())>1.5 || abs(sublead_p4.Eta())>1.5);
    TLorentzVector dipho_p4 = lead_p4 + sublead_p4;

    l.FillHist(Form("MetTag_leadGammaPt_%s",label.c_str()),    met_cat, lead_p4.Pt(), evweight);
    l.FillHist(Form("MetTag_subleadGammaPt_%s",label.c_str()), met_cat, sublead_p4.Pt(), evweight);
    l.FillHist(Form("MetTag_diphomva_%s",label.c_str()),       met_cat, bdtoutput, evweight);
    l.FillHist(Form("MetTag_uncorrmet_%s",label.c_str()),      met_cat, l.met_pfmet, evweight);
    l.FillHist(Form("MetTag_uncorrmetPhi_%s",label.c_str()),   met_cat, l.met_phi_pfmet, evweight);
    l.FillHist(Form("MetTag_corrmet_%s",label.c_str()),        met_cat, corrMet,    evweight);
    l.FillHist(Form("MetTag_corrmetPhi_%s",label.c_str()),     met_cat, corrMetPhi, evweight);
    bool cleanmet = l.METCleaning2012B(lead_p4, sublead_p4, myMet);
    if( cleanmet ) {
        l.FillHist(Form("MetTag_corrcleanedmet_%s",label.c_str()),        met_cat, corrMet,    evweight);
        l.FillHist(Form("MetTag_corrcleanedmetPhi_%s",label.c_str()),     met_cat, corrMetPhi, evweight);
    }


    if(corrMet>70){
        l.FillHist(Form("MetTag_dPhiLead_%s",label.c_str()),       met_cat, fabs(myMet.DeltaPhi(lead_p4)), evweight);
        l.FillHist(Form("MetTag_dPhiSub_%s",label.c_str()),        met_cat, fabs(myMet.DeltaPhi(sublead_p4)), evweight);
        float maxdphi = (fabs(myMet.DeltaPhi(lead_p4))>fabs(myMet.DeltaPhi(sublead_p4))) ? myMet.DeltaPhi(lead_p4) : myMet.DeltaPhi(sublead_p4);
        l.FillHist(Form("MetTag_dPhiMax_%s",label.c_str()),        met_cat, fabs(maxdphi), evweight);
        l.FillHist(Form("MetTag_dPhiDipho_%s",label.c_str()),      met_cat, fabs(myMet.DeltaPhi(dipho_p4)), evweight);
    }

    TLorentzVector* leadjt;
    float maxpt=0;
    for(int ijet=0; ijet<l.jet_algoPF1_n; ijet++){
        TLorentzVector* thisjt = (TLorentzVector*) l.jet_algoPF1_p4->At(ijet);
        if(thisjt->Pt() > maxpt){
            if(thisjt->DeltaR(lead_p4) > 0.5){
                if(thisjt->DeltaR(sublead_p4) > 0.5){
                    maxpt=thisjt->Pt();
                    leadjt=thisjt;
                }
            }
        }
    }

    if(maxpt > 50){
        l.FillHist(Form("MetTag_leadJetPt_%s",label.c_str()),      met_cat, leadjt->Pt(), evweight);
        l.FillHist(Form("MetTag_dPhiJet_%s",label.c_str()),        met_cat, fabs(myMet.DeltaPhi(*leadjt)), evweight);
        l.FillHist2D(Form("MetTag_JetPt_dPhi_%s",label.c_str()),   met_cat, leadjt->Pt(), fabs(myMet.DeltaPhi(*leadjt)), evweight);
    }

}

void PhotonAnalysis::ZWithFakeGammaCS(LoopAll& l, float* smeared_pho_energy){
    int elVtx=-1;
    bool tag = false;

    if(l.met_pfmet > 30) return;

    int el_ind=l.ElectronSelectionMVA2012();
    if(el_ind==-1) return;
    elVtx=l.FindElectronVertex(el_ind);

    int selphoind=-1;
    int selphopt=0;
    TLorentzVector* selel_p4= (TLorentzVector*) l.el_std_p4->At(el_ind);
    for(int ipho=0; ipho<l.pho_n; ipho++){
        float phoEta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[ipho]))->Eta());
        if( phoEta > 2.5 || ( phoEta > 1.4442 && phoEta < 1.566 )){
            continue;
        }
        TLorentzVector thispho = l.get_pho_p4(ipho,elVtx,smeared_pho_energy);
        if(selel_p4->DeltaR(thispho)<0.2) continue;
        float phoPt = thispho.Pt();
        if(phoPt<25) continue;
	    if (!l.PhotonMITPreSelection2011(ipho, elVtx,  smeared_pho_energy)) continue;
        float phoMVA= l.photonIDMVA(ipho, elVtx, thispho, bdtTrainingType.c_str() );
        if(phoMVA<-0.1) continue;
        if(phoPt>selphopt){
            selphopt=phoPt;
            selphoind=ipho;
            tag=true;
        }
    }

    if(!tag) return;
    TLorentzVector selpho_p4=l.get_pho_p4(selphoind,elVtx,smeared_pho_energy);

    TLorentzVector eg_p4 = selpho_p4+*selel_p4;
    int cat= (fabs(selpho_p4.Eta())<1.5 && fabs(selel_p4->Eta())<1.5);

    if(selphopt>40){
        l.FillHist("Meg_fakelead",cat,eg_p4.M());
        l.FillHist("r9_fakelead",cat,l.pho_r9[selphoind]);
    } else {
        l.FillHist("Meg_fakesub",cat,eg_p4.M());
        l.FillHist("r9_fakesub",cat,l.pho_r9[selphoind]);
    }

    return;
}


bool PhotonAnalysis::ElectronStudies2012B(LoopAll& l, float* smeared_pho_energy, bool mvaselection, float phoidMvaCut, float eventweight, float myweight, int jentry) {
    int diphotonVHlep_id=-1;
    int elVtx=-1;
    bool tag = false;

    bool debuglocal=false;

    TLorentzVector* el_tag;
    TLorentzVector* el_sc;

    float leadptcut=30;
    float subleadptcut=20;
    float elptcut=10;

    int elInd=l.ElectronSelectionMVA2012(elptcut);
    if(elInd==-1) return tag;

    if(debuglocal)  std::cout<<"ElectronStudies2012B electron"<<elInd<<std::endl;

    el_tag = (TLorentzVector*) l.el_std_p4->At(elInd);
    if(debuglocal)  std::cout<<"ElectronStudies2012B got p4"<<std::endl;
    el_sc = (TLorentzVector*) l.el_std_sc->At(elInd);
    if(debuglocal)  std::cout<<"ElectronStudies2012B got sc p4"<<std::endl;
    elVtx=l.FindElectronVertex(elInd);

    float drtoveto = 0.5;
    std::vector<bool> veto_indices;
    veto_indices.clear();
    l.PhotonsToVeto(el_sc, drtoveto, veto_indices, true);

    if(mvaselection) {
        diphotonVHlep_id = l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVHlepCut,subleadEtVHlepCut,phoidMvaCut,
            applyPtoverM, &smeared_pho_energy[0], true, true,-100, -1, false, veto_indices);
    } else {
        diphotonVHlep_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHlepCut,subleadEtVHlepCut, 4,
            applyPtoverM, &smeared_pho_energy[0], true, -1, veto_indices);
    }

    if(diphotonVHlep_id==-1) return tag;

    if(debuglocal)  std::cout<<"ElectronStudies2012B diphoton "<<diphotonVHlep_id<<std::endl;

    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);

    if(debuglocal){
        std::cout<<"lead pt "<<lead_p4.Pt()<<std::endl;
        std::cout<<"sublead pt "<<sublead_p4.Pt()<<std::endl;
        std::cout<<"el pt "<<el_tag->Pt()<<std::endl;
    }
    // need to check again for d0 and dZ (couldn't before because we didn't have the vertex)
    myEl_passelcuts = l.ElectronMVACuts(elInd, elVtx);
    myEl_ElePho     = l.ElectronPhotonCuts2012B(lead_p4, sublead_p4, *el_tag);
    if(myEl_passelcuts){
        if(myEl_ElePho){
            tag=true;
        }
    }

    if(debuglocal){
        std::cout<<"ElectronStudies2012B "<<tag<<std::endl;
    }

    int el_cat=-1;
    if(fabs(el_sc->Eta())<1.5 && fabs(lead_p4.Eta())<1.5 && fabs(sublead_p4.Eta())<1.5) el_cat=0;
    else el_cat=1;

    if(debuglocal)  std::cout<<"ElectronStudies2012B el_cat "<<el_cat<<std::endl;

    TLorentzVector gg_p4 = lead_p4 + sublead_p4;
    if(debuglocal){
        std::cout<<"el_tag->Pt(),sublead_p4.Pt(),lead_p4.Pt()/gg_p4.M() "<<el_tag->Pt()<<","<<sublead_p4.Pt()<<","<<lead_p4.Pt()/gg_p4.M()<<std::endl;
        std::cout<<"pt cuts "<<tag<<std::endl;
    }

    float thiseta = fabs(el_sc->Eta());
    float thispt  = fabs(el_tag->Pt());

    double Aeff=0.;
    if(thiseta<1.0)                   Aeff=0.10;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.12;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.085;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.11;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.12;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.12;
    if(thiseta>=2.4)                  Aeff=0.13;

    float thisiso=l.el_std_pfiso_charged[elInd]+std::max(l.el_std_pfiso_neutral[elInd]+l.el_std_pfiso_photon[elInd]-l.rho*Aeff,0.);
    TLorentzVector elpho1 = *el_tag + lead_p4;
    TLorentzVector elpho2 = *el_tag + sublead_p4;

    int cur_type = l.itype[l.current];
    if( cur_type<0 && (abs(l.gh_vh_pdgid)==24 || abs(l.gh_vh_pdgid)==26) ) {
        myEl_leptonSig = (abs(l.gh_vh1_pdgid)>10 && abs(l.gh_vh1_pdgid)<19);
    } else {
        myEl_leptonSig = 1;
    }

    float overE_overP=fabs((1/l.el_std_pin[elInd])-(1/(l.el_std_pin[elInd]*l.el_std_eopin[elInd])));
    myEl_presellead =   HLTPhotonPreselection(l, &lead_p4, elVtx);
    int matchlead=-1;
    myEl_matchellead=   PhotonMatchElectron(l, &lead_p4, matchlead);
    myEl_preselsub  =   HLTPhotonPreselection(l, &sublead_p4, elVtx);
    int matchsub=-1;
    myEl_matchelsub =   PhotonMatchElectron(l, &sublead_p4, matchsub);
    myEl_ptlead     =   lead_p4.Pt();
    myEl_ptsub      =   sublead_p4.Pt();
    myEl_elvetolead =   l.pho_isconv[l.dipho_leadind[diphotonVHlep_id]];
    myEl_elvetosub  =   l.pho_isconv[l.dipho_subleadind[diphotonVHlep_id]];

    int lead_cat = abs(lead_p4.Eta())>1.5;
    int sublead_cat = abs(sublead_p4.Eta())>1.5;
    l.FillHist("PhoMatchEl",lead_cat,   myEl_matchellead);
    l.FillHist("PhoMatchEl",sublead_cat,myEl_matchelsub);

    if(myEl_matchellead==1){
        if(cur_type==0) std::cout<<"badphoton event,lumis,run:  "<<l.event<<","<<l.lumis<<","<<l.run<<std::endl;
        l.FillHist("MatchedPhoPt",      lead_cat, myEl_ptlead);
        l.FillHist("MatchedPhoElVeto",  lead_cat, myEl_elvetolead);
        l.FillHist("MatchedPhoElTrk",   lead_cat, l.el_std_tkdrv[matchlead]);
        l.FillHist("MatchedPhoElEcal",  lead_cat, l.el_std_ecaldrv[matchlead]);
        l.FillHist("MatchedPhoGenPho",  lead_cat, l.pho_genmatched[l.dipho_leadind[diphotonVHlep_id]]);
        l.FillHist("MatchedSCIndex",    lead_cat, (int)(l.el_std_scind[matchlead]==l.pho_scind[l.dipho_leadind[diphotonVHlep_id]]));
    }

    if(myEl_matchelsub==1){
        if(cur_type==0) std::cout<<"badphoton event,lumis,run:  "<<l.event<<","<<l.lumis<<","<<l.run<<std::endl;
        l.FillHist("MatchedPhoPt",      sublead_cat, myEl_ptsub);
        l.FillHist("MatchedPhoElVeto",  sublead_cat, myEl_elvetosub);
        l.FillHist("MatchedPhoElTrk",   sublead_cat, l.el_std_tkdrv[matchsub]);
        l.FillHist("MatchedPhoElEcal",  sublead_cat, l.el_std_ecaldrv[matchsub]);
        l.FillHist("MatchedPhoGenPho",  sublead_cat, l.pho_genmatched[l.dipho_subleadind[diphotonVHlep_id]]);
        l.FillHist("MatchedSCIndex",    sublead_cat, (int)(l.el_std_scind[matchsub]==l.pho_scind[l.dipho_subleadind[diphotonVHlep_id]]));
    }

    if(debuglocal)  std::cout<<"ElectronStudies2012B filled pho variables "<<std::endl;

    myEl_elpt       =   thispt;
    myEl_oEsuboP    =   overE_overP;
    myEl_D0         =   l.el_std_D0Vtx[elInd][elVtx];
    myEl_DZ         =   l.el_std_DZVtx[elInd][elVtx];
    myEl_conv       =   l.el_std_conv[elInd];
    myEl_mishit     =   l.el_std_hp_expin[elInd];
    myEl_reliso     =   thisiso/thispt;
    myEl_iso        =   thisiso;
    myEl_detain     =   l.el_std_detain[elInd];
    myEl_dphiin     =   l.el_std_dphiin[elInd];
    myEl_sieie      =   l.el_std_sieie[elInd];
    myEl_hoe        =   l.el_std_hoe[elInd];

    if(debuglocal)  std::cout<<"ElectronStudies2012B filled pho+ele variables "<<std::endl;

    myEl_drlead     =   lead_p4.DeltaR(*el_tag);
    myEl_drsub      =   sublead_p4.DeltaR(*el_tag);
    myEl_melead     =   elpho1.M();
    myEl_meleadveto15 = abs(myEl_melead-91.2)<15;
    myEl_meleadveto10 = abs(myEl_melead-91.2)<10;
    myEl_mesub      =   elpho2.M();
    myEl_mesubveto10=   abs(myEl_mesub -91.2)<10;
    myEl_mesubveto5 =   abs(myEl_mesub -91.2)<5;
    //myEl_mvaTrig    =   l.el_std_mva_trig[elInd];
    myEl_mvaNonTrig =   l.el_std_mva_nontrig[elInd];

    //std::cout<<"myEl_ElePho "<<myEl_ElePho<<std::endl;

    myEl_ptgg       =   gg_p4.Pt();
    myEl_mgg        =   gg_p4.M();
    myEl_ptleadom   =   myEl_ptlead/myEl_mgg;
    myEl_ptsubom    =   myEl_ptsub/myEl_mgg;
    myEl_phomaxeta  =   std::max(lead_p4.Eta(),sublead_p4.Eta());
    myEl_sumpt3     =   lead_p4.Pt()+sublead_p4.Pt()+el_tag->Pt();
    myEl_dRtklead   =   l.pho_drtotk_25_99[l.dipho_leadind[diphotonVHlep_id]];
    myEl_dRtksub    =   l.pho_drtotk_25_99[l.dipho_subleadind[diphotonVHlep_id]];

    if( debuglocal ) {
        if(myEl_dRtklead<1.0 || myEl_dRtksub<1.0 || myEl_matchellead==1 || myEl_matchelsub==1){
            std::cout<<"myEl_dRtklead,myEl_dRtksub,myEl_matchellead,myEl_matchelsub "<<
                    myEl_dRtklead<<","<<
                    myEl_dRtksub<<","<<
                    myEl_matchellead<<","<<
                    myEl_matchelsub<<std::endl;
        }
    }

    myEl_MVAlead    = l.photonIDMVA(l.dipho_leadind[diphotonVHlep_id], elVtx,
                                    lead_p4,bdtTrainingType.c_str());
    myEl_MVAsub     = l.photonIDMVA(l.dipho_subleadind[diphotonVHlep_id], elVtx,
                                        sublead_p4,bdtTrainingType.c_str());

    if(debuglocal && myEl_ElePho){
        std::cout<<"myEl_mvaNonTrig,myEl_MVAlead,myEl_MVAsub "
            <<myEl_mvaNonTrig<<","<<myEl_MVAlead<<","<<myEl_MVAsub<<std::endl;
    }

    if( debuglocal ) std::cout<<"test01"<<std::endl;

    // Mass Resolution of the Event
    massResolutionCalculator->Setup(l,&photonInfoCollection[l.dipho_leadind[diphotonVHlep_id]],&photonInfoCollection[l.dipho_subleadind[diphotonVHlep_id]],elVtx,massResoPars,nR9Categories,nEtaCategories,beamspotSigma,true);
//    float vtx_mva  = l.vtx_std_evt_mva->at(diphotonVHlep_id);
    //for(int rankind=0; rankind<l.vtx_std_ranked_list->size(); rankind++){
    //    if(l.vtx_std_ranked_list[rankind]==elVtx){
    //        l.vtx_std_ranked_list->erase(l.vtx_std_ranked_list->begin()+rankind);
    //    }
    //}
    //l.vtx_std_ranked_list->insert(l.vtx_std_ranked_list->begin(),elVtx);
    float vtx_mva = -0.9;// vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back()  );
    float sigmaMrv = massResolutionCalculator->massResolutionEonly();
    float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
    float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
    // easy to calculate vertex probability from vtx mva output
    float vtxProb   = 1.-0.49*(vtx_mva+1.0); /// should better use this: vtxAna_.setPairID(diphoton_id); vtxAna_.vertexProbability(vtx_mva); PM
    if( debuglocal ) std::cout<<"test02"<<std::endl;

	float diphobdt_output = l.diphotonMVA(-1,l.dipho_leadind[diphotonVHlep_id], l.dipho_subleadind[diphotonVHlep_id], elVtx,
                                          vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,
                                          bdtTrainingPhilosophy.c_str(), bdtTrainingType.c_str(),
                                          myEl_MVAlead,myEl_MVAsub);

    myEl_diphomva    =   diphobdt_output;

    std::vector<std::vector<bool> > ph_passcut;
    if( debuglocal ) std::cout<<"test03"<<std::endl;
    myEl_CiClead    =   l.PhotonCiCPFSelectionLevel(l.dipho_leadind[diphotonVHlep_id], elVtx, ph_passcut, 4, 0, smeared_pho_energy );
    if( debuglocal ) std::cout<<"test04"<<std::endl;
    myEl_CiCsub     =   l.PhotonCiCPFSelectionLevel(l.dipho_subleadind[diphotonVHlep_id], elVtx, ph_passcut, 4, 1, smeared_pho_energy );
    if( debuglocal ) std::cout<<"test05"<<std::endl;
    myEl_MET        =   l.met_pfmet;
    myEl_METphi     =   l.met_phi_pfmet;

    int elcat = (thiseta>1.5 || fabs(lead_p4.Eta())>1.5 || fabs(sublead_p4.Eta())>1.5);
    if(myEl_diphomva>0.05 && (myEl_meleadveto15 || myEl_mesubveto10)) {
        l.FillHist("NearZ_leadr9",    elcat,   (float)l.pho_r9[l.dipho_leadind[diphotonVHlep_id]]);
        if(l.pho_r9[l.dipho_leadind[diphotonVHlep_id]]>0.94){
            l.FillHist("M_elead_hir9",    elcat,   myEl_melead);
        } else {
            l.FillHist("M_elead_lor9",    elcat,   myEl_melead);
        }
        l.FillHist("NearZ_subleadr9",    elcat,   (float)l.pho_r9[l.dipho_subleadind[diphotonVHlep_id]]);
        if(l.pho_r9[l.dipho_subleadind[diphotonVHlep_id]]>0.94){
            l.FillHist("M_esublead_hir9",    elcat,   myEl_mesub);
        } else {
            l.FillHist("M_esublead_lor9",    elcat,   myEl_mesub);
        }
    }
    myEl_category=(int)(fabs(lead_p4.Eta())>1.5)+2*(int)(fabs(sublead_p4.Eta())>1.5) +4*(int)(thiseta>1.5);
    if( debuglocal ) std::cout<<"test08"<<std::endl;
    tag = l.ApplyCutsFill(elcat,10, eventweight, myweight);
    if(tag&&cur_type==0) std::cout<<"selected electron tag event,lumis,run:  "<<l.event<<","<<l.lumis<<","<<l.run<<std::endl;
    if( debuglocal ) std::cout<<"test09"<<std::endl;

    //plots for comparing selection

    if(diphotonVHlep_id!=-1){
        l.FillHist("ElectronTag_MethodA_MethodB_samevtx", 0,(float)(l.dipho_vtxind[diphotonVHlep_id]==elVtx));
    }

    return tag;
}


bool PhotonAnalysis::ElectronTagStudies2012(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight, int jentry){
    bool tag = false;

    bool debuglocal=false;

    TLorentzVector* el_tag;
    TLorentzVector* el_sc;
    int elVtx = 0;

    TLorentzVector lead_p4;
    int leadpho_ind=-1;

    TLorentzVector sublead_p4;
    int subleadpho_ind=-1;

    float leadptcut=30;
    float subleadptcut=20;
    float elptcut=10;

    int elInd=l.ElectronSelectionMVA2012(elptcut);
    if(elInd!=-1) {
        el_tag = (TLorentzVector*) l.el_std_p4->At(elInd);
        if(el_tag->Pt()>elptcut){
            el_sc = (TLorentzVector*) l.el_std_sc->At(elInd);
            elVtx=l.FindElectronVertex(elInd);

            float drtoveto = 0.2;
            std::vector<bool> veto_indices;
            veto_indices.clear();
            l.PhotonsToVeto(el_sc, drtoveto, veto_indices, true);

            diphotonVHlep_id = l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadptcut,subleadptcut,-0.2,
                    applyPtoverM, &smeared_pho_energy[0], true, true,-100, -1, false, veto_indices);
            //    diphotonVHlep_id = l.DiphotonCiCSelection( l.phoLOOSE, l.phoLOOSE, leadEtVHlepCut,subleadEtVHlepCut, 4,
            //        applyPtoverM, &smeared_pho_energy[0], true, elVtx, veto_indices);

            if(diphotonVHlep_id!=-1){
                lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], elVtx, &smeared_pho_energy[0]);
                sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], elVtx, &smeared_pho_energy[0]);
                tag=true;
            }
        }
    }

    if(!tag) return tag;

    if( debuglocal ) {
        std::cout<<"l.ElectronPhotonCuts(lead_p4, sublead_p4, *el_tag)"<<l.ElectronPhotonCuts(lead_p4, sublead_p4, *el_tag)<<std::endl;
    }
    // require tighter dr to electron track
    //if(l.pho_drtotk_25_99[l.dipho_leadind[diphotonVHlep_id]] < 1 || l.pho_drtotk_25_99[l.dipho_subleadind[diphotonVHlep_id]] < 1) return tag;




    float thiseta = fabs(el_sc->Eta());
    float thispt  = fabs(el_tag->Pt());

    double Aeff=0.;
    if(thiseta<1.0)                   Aeff=0.10;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.12;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.085;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.11;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.12;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.12;
    if(thiseta>=2.4)                  Aeff=0.13;

    float thisiso=l.el_std_pfiso_charged[elInd]+std::max(l.el_std_pfiso_neutral[elInd]+l.el_std_pfiso_photon[elInd]-l.rho*Aeff,0.);
    TLorentzVector elpho1 = *el_tag + lead_p4;
    TLorentzVector elpho2 = *el_tag + sublead_p4;

    int cur_type = l.itype[l.current];
    if( cur_type<0 && (abs(l.gh_vh_pdgid)==24 || abs(l.gh_vh_pdgid)==26) ) {
        myEl_leptonSig = (abs(l.gh_vh1_pdgid)>10 && abs(l.gh_vh1_pdgid)<19);
    } else {
        myEl_leptonSig = 1;
    }

    float overE_overP=fabs((1/l.el_std_pin[elInd])-(1/(l.el_std_pin[elInd]*l.el_std_eopin[elInd])));
    myEl_presellead =   HLTPhotonPreselection(l, &lead_p4, elVtx);
    int matchlead=-1;
    myEl_matchellead=   PhotonMatchElectron(l, &lead_p4, matchlead);
    myEl_preselsub  =   HLTPhotonPreselection(l, &sublead_p4, elVtx);
    int matchsub=-1;
    myEl_matchelsub =   PhotonMatchElectron(l, &sublead_p4, matchsub);
    myEl_ptlead     =   lead_p4.Pt();
    myEl_ptsub      =   sublead_p4.Pt();
    myEl_elvetolead =   l.pho_isconv[leadpho_ind];
    myEl_elvetosub  =   l.pho_isconv[subleadpho_ind];

    int lead_cat = abs(lead_p4.Eta())>1.5;
    int sublead_cat = abs(sublead_p4.Eta())>1.5;
    l.FillHist("PhoMatchEl",lead_cat,   myEl_matchellead);
    l.FillHist("PhoMatchEl",sublead_cat,myEl_matchelsub);

    if(myEl_matchellead==1){
        if(cur_type==0) std::cout<<"badphoton event,lumis,run:  "<<l.event<<","<<l.lumis<<","<<l.run<<std::endl;
        l.FillHist("MatchedPhoPt",      lead_cat, myEl_ptlead);
        l.FillHist("MatchedPhoElVeto",  lead_cat, myEl_elvetolead);
        l.FillHist("MatchedPhoElTrk",   lead_cat, l.el_std_tkdrv[matchlead]);
        l.FillHist("MatchedPhoElEcal",  lead_cat, l.el_std_ecaldrv[matchlead]);
        l.FillHist("MatchedPhoGenPho",  lead_cat, l.pho_genmatched[leadpho_ind]);
        l.FillHist("MatchedSCIndex",    lead_cat, (int)(l.el_std_scind[matchlead]==l.pho_scind[leadpho_ind]));
    }

    if(myEl_matchelsub==1){
        if(cur_type==0) std::cout<<"badphoton event,lumis,run:  "<<l.event<<","<<l.lumis<<","<<l.run<<std::endl;
        l.FillHist("MatchedPhoPt",      sublead_cat, myEl_ptsub);
        l.FillHist("MatchedPhoElVeto",  sublead_cat, myEl_elvetosub);
        l.FillHist("MatchedPhoElTrk",   sublead_cat, l.el_std_tkdrv[matchsub]);
        l.FillHist("MatchedPhoElEcal",  sublead_cat, l.el_std_ecaldrv[matchsub]);
        l.FillHist("MatchedPhoGenPho",  sublead_cat, l.pho_genmatched[subleadpho_ind]);
        l.FillHist("MatchedSCIndex",    sublead_cat, (int)(l.el_std_scind[matchsub]==l.pho_scind[subleadpho_ind]));
    }

    myEl_elpt       =   thispt;
    myEl_oEsuboP    =   overE_overP;
    myEl_D0         =   l.el_std_D0Vtx[elInd][elVtx];
    myEl_DZ         =   l.el_std_DZVtx[elInd][elVtx];
    myEl_conv       =   l.el_std_conv[elInd];
    myEl_mishit     =   l.el_std_hp_expin[elInd];
    myEl_reliso     =   thisiso/thispt;
    myEl_iso        =   thisiso;
    myEl_detain     =   l.el_std_detain[elInd];
    myEl_dphiin     =   l.el_std_dphiin[elInd];
    myEl_sieie      =   l.el_std_sieie[elInd];
    myEl_hoe        =   l.el_std_hoe[elInd];

    myEl_drlead     =   lead_p4.DeltaR(*el_tag);
    myEl_drsub      =   sublead_p4.DeltaR(*el_tag);
    myEl_melead     =   elpho1.M();
    myEl_meleadveto15 = abs(myEl_melead-91.2)<15;
    myEl_meleadveto10 = abs(myEl_melead-91.2)<10;
    myEl_mesub      =   elpho2.M();
    myEl_mesubveto10=   abs(myEl_mesub -91.2)<10;
    myEl_mesubveto5 =   abs(myEl_mesub -91.2)<5;
    //myEl_mvaTrig    =   l.el_std_mva_trig[elInd];
    myEl_mvaNonTrig =   l.el_std_mva_nontrig[elInd];

    myEl_ElePho     =   l.ElectronPhotonCuts2012B(lead_p4, sublead_p4, *el_tag);
    //std::cout<<"myEl_ElePho "<<myEl_ElePho<<std::endl;

    TLorentzVector gg_p4 = lead_p4+sublead_p4;
    myEl_ptgg       =   gg_p4.Pt();
    myEl_mgg        =   gg_p4.M();
    myEl_ptleadom   =   myEl_ptlead/myEl_mgg;
    myEl_ptsubom    =   myEl_ptsub/myEl_mgg;
    myEl_phomaxeta  =   std::max(lead_p4.Eta(),sublead_p4.Eta());
    myEl_sumpt3     =   lead_p4.Pt()+sublead_p4.Pt()+el_tag->Pt();
    myEl_dRtklead   =   l.pho_drtotk_25_99[leadpho_ind];
    myEl_dRtksub    =   l.pho_drtotk_25_99[subleadpho_ind];

    if( debuglocal ) {
        if(myEl_dRtklead<1.0 || myEl_dRtksub<1.0 || myEl_matchellead==1 || myEl_matchelsub==1){
            std::cout<<"myEl_dRtklead,myEl_dRtksub,myEl_matchellead,myEl_matchelsub "<<
                    myEl_dRtklead<<","<<
                    myEl_dRtksub<<","<<
                    myEl_matchellead<<","<<
                    myEl_matchelsub<<std::endl;
        }
    }

    myEl_MVAlead    = l.photonIDMVA(leadpho_ind, elVtx,
                                    lead_p4,bdtTrainingType.c_str());
    myEl_MVAsub     = l.photonIDMVA(subleadpho_ind, elVtx,
                                        sublead_p4,bdtTrainingType.c_str());

    if(debuglocal && myEl_ElePho){
        std::cout<<"myEl_mvaNonTrig,myEl_MVAlead,myEl_MVAsub "
            <<myEl_mvaNonTrig<<","<<myEl_MVAlead<<","<<myEl_MVAsub<<std::endl;
    }

    if( debuglocal ) std::cout<<"test01"<<std::endl;

    // Mass Resolution of the Event
    massResolutionCalculator->Setup(l,&photonInfoCollection[leadpho_ind],&photonInfoCollection[subleadpho_ind],elVtx,massResoPars,nR9Categories,nEtaCategories,beamspotSigma,true);
//    float vtx_mva  = l.vtx_std_evt_mva->at(diphotonVHlep_id);
    //for(int rankind=0; rankind<l.vtx_std_ranked_list->size(); rankind++){
    //    if(l.vtx_std_ranked_list[rankind]==elVtx){
    //        l.vtx_std_ranked_list->erase(l.vtx_std_ranked_list->begin()+rankind);
    //    }
    //}
    //l.vtx_std_ranked_list->insert(l.vtx_std_ranked_list->begin(),elVtx);
    float vtx_mva = -0.9;// vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back()  );
    float sigmaMrv = massResolutionCalculator->massResolutionEonly();
    float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
    float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
    // easy to calculate vertex probability from vtx mva output
    float vtxProb   = 1.-0.49*(vtx_mva+1.0); /// should better use this: vtxAna_.setPairID(diphoton_id); vtxAna_.vertexProbability(vtx_mva); PM
    if( debuglocal ) std::cout<<"test02"<<std::endl;

	float diphobdt_output = l.diphotonMVA(-1,leadpho_ind, subleadpho_ind, elVtx,
                                          vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,
                                          bdtTrainingPhilosophy.c_str(), bdtTrainingType.c_str(),
                                          myEl_MVAlead,myEl_MVAsub);

    myEl_diphomva    =   diphobdt_output;

    std::vector<std::vector<bool> > ph_passcut;
    if( debuglocal ) std::cout<<"test03"<<std::endl;
    myEl_CiClead    =   l.PhotonCiCPFSelectionLevel(leadpho_ind, elVtx, ph_passcut, 4, 0, smeared_pho_energy );
    if( debuglocal ) std::cout<<"test04"<<std::endl;
    myEl_CiCsub     =   l.PhotonCiCPFSelectionLevel(subleadpho_ind, elVtx, ph_passcut, 4, 1, smeared_pho_energy );
    if( debuglocal ) std::cout<<"test05"<<std::endl;
    myEl_MET        =   l.met_pfmet;
    myEl_METphi     =   l.met_phi_pfmet;


    int elcat = (thiseta>1.5 || fabs(lead_p4.Eta())>1.5 || fabs(sublead_p4.Eta())>1.5);
    myEl_category=(int)(fabs(lead_p4.Eta())>1.5)+2*(int)(fabs(sublead_p4.Eta())>1.5) +4*(int)(thiseta>1.5);
    if( debuglocal ) std::cout<<"test08"<<std::endl;
    tag = l.ApplyCutsFill(elcat,10, eventweight, myweight);
    if(tag&&cur_type==0) std::cout<<"selected electron tag event,lumis,run:  "<<l.event<<","<<l.lumis<<","<<l.run<<std::endl;
    if( debuglocal ) std::cout<<"test09"<<std::endl;
    //if(jentry>-1) {
    //    if(cur_type!=0){
    //        std::map<std::string,int> genbranchnames;
    //        const char* branchchar[] = {"gp_n",
    //            "gp_mother",
    //            "gp_status",
    //            "gp_pdgid",
    //            "gp_p4"};
    //        std::vector<std::string> branchstrings(branchchar,branchchar + 5);;
    //
    //        for (int istr=0; istr<branchstrings.size(); istr++){
    //            genbranchnames[branchstrings[istr]]=1;
    //        }
    //
    //        std::set<TBranch *> genbranches;
    //        genbranches.clear();
    //        l.GetBranches(genbranchnames, genbranches);
    //        l.SetBranchAddresses(genbranchnames);
    //        l.GetEntry(genbranches, jentry);
    //    }
    //
    //    int el_gen_ind = -1;
    //    if( debuglocal ) std::cout<<"about to gen match el"<<std::endl;
    //    if(cur_type!=0){
    //        el_gen_ind=GenMatch(l, el_tag);
    //    }
    //    l.FillTree("elGenIndex",(float)el_gen_ind);
    //    l.FillTree("elEta",     (float)el_tag->Eta());
    //    l.FillTree("elPhi",     (float)el_tag->Phi());
    //    l.FillTree("elPt",      (float)el_tag->Pt());
    //    if(el_gen_ind!=-1){
    //        TLorentzVector* genel_p4=(TLorentzVector*) l.gp_p4->At(el_gen_ind);
    //        l.FillTree("elGenEta",  (float)genel_p4->Eta());
    //        l.FillTree("elGenPhi",  (float)genel_p4->Phi());
    //        l.FillTree("elGenPt",   (float)genel_p4->Pt());
    //        l.FillTree("elGenPDGID",(float)l.gp_pdgid[el_gen_ind]);
    //        float dr_min = el_tag->DeltaR(*genel_p4);
    //        float dpt_min = el_tag->Pt() - genel_p4->Pt();
    //        l.FillTree("elGenDR",   (float)dr_min);
    //        l.FillTree("elGenDPT",  (float)dpt_min);
    //    }
    //    l.FillTree("elmva", (float)myEl_mvaNonTrig);

    //    int pholead_gen_ind = -1;
    //    if(cur_type!=0){
    //        pholead_gen_ind=GenMatch(l, &lead_p4);
    //    }
    //    l.FillTree("leadGenIndex",(float)pholead_gen_ind);
    //    l.FillTree("leadEta",     (float)lead_p4.Eta());
    //    l.FillTree("leadPhi",     (float)lead_p4.Phi());
    //    l.FillTree("leadPt",      (float)lead_p4.Pt());
    //    l.FillTree("leadr9",      (float)l.pho_r9[leadpho_ind]);
    //    if(pholead_gen_ind!=-1){
    //        TLorentzVector* genlead_p4=(TLorentzVector*) l.gp_p4->At(pholead_gen_ind);
    //        l.FillTree("leadGenEta",  (float)genlead_p4->Eta());
    //        l.FillTree("leadGenPhi",  (float)genlead_p4->Phi());
    //        l.FillTree("leadGenPt",   (float)genlead_p4->Pt());
    //        l.FillTree("leadGenPDGID",(float)l.gp_pdgid[pholead_gen_ind]);
    //        float dr_min  = lead_p4.DeltaR(*genlead_p4);
    //        float dpt_min = lead_p4.Pt() - genlead_p4->Pt();
    //        l.FillTree("leadGenDR",   (float)dr_min);
    //        l.FillTree("leadGenDPT",  (float)dpt_min);
    //    }
    //    l.FillTree("lead_eleveto", (float)myEl_elvetolead);
    //    l.FillTree("sub_eleveto", (float)myEl_elvetosub);
    //    l.FillTree("leadCIC", (float)myEl_CiClead);
    //    l.FillTree("subCIC", (float)myEl_CiCsub);
    //    l.FillTree("lead_Elmatched", (float)myEl_matchellead);
    //    l.FillTree("sub_Elmatched", (float)myEl_matchelsub);
    //    l.FillTree("leadmva", (float)myEl_MVAlead);
    //    l.FillTree("submva", (float)myEl_MVAsub);

    //    int phosublead_gen_ind=-1;
    //    if(cur_type!=0){
    //        phosublead_gen_ind=GenMatch(l, &sublead_p4);
    //    }
    //    l.FillTree("subleadGenIndex",(float)phosublead_gen_ind);
    //    l.FillTree("subleadEta",     (float)sublead_p4.Eta());
    //    l.FillTree("subleadPhi",     (float)sublead_p4.Phi());
    //    l.FillTree("subleadPt",      (float)sublead_p4.Pt());
    //    l.FillTree("subleadr9",      (float)l.pho_r9[subleadpho_ind]);
    //    if(phosublead_gen_ind!=-1){
    //        TLorentzVector* gensublead_p4=(TLorentzVector*) l.gp_p4->At(phosublead_gen_ind);
    //        l.FillTree("subleadGenEta",  (float)gensublead_p4->Eta());
    //        l.FillTree("subleadGenPhi",  (float)gensublead_p4->Phi());
    //        l.FillTree("subleadGenPt",   (float)gensublead_p4->Pt());
    //        l.FillTree("subleadGenPDGID",(float)l.gp_pdgid[phosublead_gen_ind]);
    //        float dr_min  = sublead_p4.DeltaR(*gensublead_p4);
    //        float dpt_min = sublead_p4.Pt() - gensublead_p4->Pt();
    //        l.FillTree("subleadGenDR",   (float)dr_min);
    //        l.FillTree("subleadGenDPT",  (float)dpt_min);
    //    }


    //    l.FillTree("dipho_mva", (float)myEl_diphomva);
    //    l.FillTree("cur_type",  (int)cur_type);
    //
    //    l.FillTree("tagged",    (float)tag);
    //    l.FillTree("weight",    (float)eventweight);
    //}

    //plots for comparing selection

    if(diphotonVHlep_id!=-1){
        l.FillHist("ElectronTag_MethodA_MethodB_samelead",0,(float)(l.dipho_leadind[diphotonVHlep_id]==leadpho_ind));
        l.FillHist("ElectronTag_MethodA_MethodB_samesub", 0,(float)(l.dipho_subleadind[diphotonVHlep_id]==subleadpho_ind));
        l.FillHist("ElectronTag_MethodA_MethodB_samevtx", 0,(float)(l.dipho_vtxind[diphotonVHlep_id]==elVtx));
    }

    return tag;
}



int PhotonAnalysis::GenMatch(LoopAll& l, TLorentzVector* recop4){

    int reco_gen_ind=-1;
    int reco_gen_pdgid=-1000;
    TLorentzVector* igp_p4;
    TLorentzVector* genmatch_p4;
    float dr_temp=999;
    float dr_min=999;
    float dpt_temp=5;
    float dpt_min=5;
    for(int igp=0; igp<l.gp_n; igp++){
        if(l.gp_status[igp] != 1) continue;
        igp_p4 = (TLorentzVector*) l.gp_p4->At(igp);
        dr_temp = recop4->DeltaR(*igp_p4);
        dpt_temp = recop4->Pt() - igp_p4->Pt();
        if(dr_temp<dr_min && abs(dpt_temp)*1.15<abs(dpt_min)){
            dr_min=dr_temp;
            dpt_min=dpt_temp;
            reco_gen_ind=igp;
            reco_gen_pdgid=l.gp_pdgid[igp];
            genmatch_p4=(TLorentzVector*) l.gp_p4->At(igp);
        }
        if(dr_min<0.05) break;
    }

    return reco_gen_ind;
}

bool PhotonAnalysis::HLTPhotonPreselection(LoopAll& l, TLorentzVector* phop4, int photon_index){
    bool pass = true;
    // denominator definition: only matching HLT related part of MIT preselection

    // cuts to emulate the trigger
    float HILTcuts_hoe[4]     = {0.082,0.075,0.075,0.075};
    float HILTcuts_sieie[4]   = {0.014,0.014,0.034,0.034};
    float HILTcuts_ecaliso[4] = {50,4,50,4};
    float HILTcuts_hcaliso[4] = {50,4,50,4};
    float HILTcuts_trkiso[4]  = {50,4,50,4};

    // pT cut
    float phopt = phop4->Et();
    if (phopt<20) pass = false;

    // eta cut
    if( phop4->Eta() > 2.5 || ( phop4->Eta()>1.4442 && phop4->Eta()<1.566 ) ) pass = false;

    // photon category
    int PhotonEtaCategory;
    if (l.pho_isEB[photon_index]) PhotonEtaCategory=0;
    if (l.pho_isEE[photon_index]) PhotonEtaCategory=1;
    int r9_category     = (int) (l.pho_r9[photon_index] <= 0.9);
    int photon_category = r9_category + 2*PhotonEtaCategory;

    // isolation and id cuts
    float val_hoe        = l.pho_hoe[photon_index];
    float val_sieie      = l.pho_sieie[photon_index];
    float val_ecaliso    = l.pho_ecalsumetconedr03[photon_index] - 0.012*phop4->Et();
    float val_hcaliso    = l.pho_hcalsumetconedr03[photon_index] - 0.005*phop4->Et();
    float val_trkiso     = l.pho_trksumpthollowconedr03[photon_index] - 0.002*phop4->Et();
    if (val_hoe     >= HILTcuts_hoe[photon_category]     ) pass = false;
    if (val_sieie   >= HILTcuts_sieie[photon_category]   ) pass = false;
    if (val_ecaliso >= HILTcuts_ecaliso[photon_category] ) pass = false;
    if (val_hcaliso >= HILTcuts_hcaliso[photon_category] ) pass = false;
    if (val_trkiso  >= HILTcuts_trkiso[photon_category]  ) pass = false;

    int val_pho_isconv = l.pho_isconv[photon_index];
    if ( !val_pho_isconv ) pass = false;

    return pass;
}


bool PhotonAnalysis::PhotonMatchElectron(LoopAll& l, TLorentzVector* pho_p4){
    int dummy=-1;

    return PhotonMatchElectron(l, pho_p4, dummy);
}

bool PhotonAnalysis::PhotonMatchElectron(LoopAll& l, TLorentzVector* pho_p4, int& el_match_ind){
    bool pass = false;

    // need to check if this photon matches an electron with 0 missing hits
    el_match_ind=-1;
    TLorentzVector* el_p4;
    for(int iel=0; iel<l.el_std_n; iel++){
        if(l.el_std_conv[iel]==1 || l.el_std_hp_expin[iel]!=0) continue;
        el_p4=(TLorentzVector*) l.el_std_sc->At(iel);
        float dr = el_p4->DeltaR(*pho_p4);
        float dpt = el_p4->Pt() - pho_p4->Pt();
        if(dr<0.2 && abs(dpt/pho_p4->Pt())<0.05) {
            el_match_ind=iel;
            pass=true;
            break;
        }
    }
    return pass;
}


bool PhotonAnalysis::MuonTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphotonVHlep_id==-1) return tag;

    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);

    int muonInd = l.MuonSelection(lead_p4, sublead_p4, l.dipho_vtxind[diphotonVHlep_id]);
    if(muonInd!=-1) tag = true;

    return tag;
}

bool PhotonAnalysis::MuonTag2012(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, ofstream& lep_sync, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphotonVHlep_id==-1) return tag;

    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);

    int muonInd = l.MuonSelection2012(lead_p4, sublead_p4, l.dipho_vtxind[diphotonVHlep_id]);
    if(muonInd!=-1) tag = true;

        TLorentzVector dipho_p4 = lead_p4+sublead_p4;
        lep_sync<<"run="<<l.run<<"\t";
        lep_sync<<"lumis"<<l.lumis<<"\t";
        lep_sync<<"event="<<(unsigned int) l.event<<"\t";
        lep_sync<<"pt1="<<lead_p4.Pt()<<"\t";
        lep_sync<<"eta1="<<lead_p4.Eta()<<"\t";
        lep_sync<<"pt2="<<sublead_p4.Pt()<<"\t";
        lep_sync<<"eta2="<<sublead_p4.Eta()<<"\t";
        lep_sync<<"ptgg="<<dipho_p4.Pt()<<"\t";
    if(tag){
        TLorentzVector* mu_p4 = (TLorentzVector*) l.mu_glo_p4->At(muonInd);
        lep_sync<<"mupt="<<mu_p4->Pt()<<"\t";
        lep_sync<<"mueta="<<mu_p4->Eta()<<"  MUTAG\n";
    } else {
        lep_sync<<"MUNOTAG\n";
    }

    return tag;
}

bool PhotonAnalysis::MuonTag2012B(LoopAll& l, int& diphotonVHlep_id, int& mu_ind, int& muVtx, int& mu_cat, float* smeared_pho_energy, ofstream& lep_sync, bool mvaselection, float phoidMvaCut, float eventweight, std::vector<float> smeared_pho_weight, bool fillHist, bool vetodipho, bool kinonly){
    bool tag = false;
    float muptcut=20.;

    mu_ind=l.MuonSelection2012B(muptcut);

    if(mu_ind!=-1) {
        TLorentzVector* mymu = (TLorentzVector*) l.mu_glo_p4->At(mu_ind);
        muVtx=l.FindMuonVertex(mu_ind);

        float drtoveto = 0.5;
        std::vector<bool> veto_indices;
        veto_indices.clear();
        l.PhotonsToVeto(mymu, drtoveto, veto_indices, false);

        if(mvaselection) {
            diphotonVHlep_id = l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVHlepCut,subleadEtVHlepCut,phoidMvaCut,
                applyPtoverM, &smeared_pho_energy[0], vetodipho, kinonly, diphobdt_output_Cut_VHLep, -1, false, veto_indices);
        } else {
            diphotonVHlep_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHlepCut,subleadEtVHlepCut, 4,
                applyPtoverM, &smeared_pho_energy[0], true, -1, veto_indices);
        }

        if(diphotonVHlep_id!=-1){
            TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
            TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
            TLorentzVector dipho_p4 = lead_p4 + sublead_p4;
            float mass = dipho_p4.M();
            std::string label("noleppho_nomva");
            if(mass>=100 && mass<180 && fillHist){
                if(smeared_pho_weight.size()!=0){
                    eventweight*=(smeared_pho_weight[l.dipho_leadind[diphotonVHlep_id]] * smeared_pho_weight[l.dipho_subleadind[diphotonVHlep_id]]);
                }
                int cur_type = l.itype[l.current];
                ControlPlotsMuonTag2012B(l, lead_p4, sublead_p4, mu_ind, 0, eventweight, label);
            }

            tag = l.MuonPhotonCuts2012B(lead_p4, sublead_p4, mymu);
            if(!tag) diphotonVHlep_id=-1;
            mu_cat=(int)(abs(lead_p4.Eta())>1.5 || abs(sublead_p4.Eta())>1.5);
            if(tag){
                lep_sync<<"run="<<l.run<<"\t";
                lep_sync<<"lumis"<<l.lumis<<"\t";
                lep_sync<<"event="<<(unsigned int) l.event<<"\t";
                lep_sync<<"pt1="<<lead_p4.Pt()<<"\t";
                lep_sync<<"eta1="<<lead_p4.Eta()<<"\t";
                lep_sync<<"pt2="<<sublead_p4.Pt()<<"\t";
                lep_sync<<"eta2="<<sublead_p4.Eta()<<"\t";
                lep_sync<<"ptgg="<<dipho_p4.Pt()<<"\t";
                lep_sync<<"mupt="<<mymu->Pt()<<"\t";
                lep_sync<<"mueta="<<mymu->Eta()<<"  MUTAG\n";
            }
        }
    }

    return tag;
}


void PhotonAnalysis::ControlPlotsMuonTag2012B(LoopAll& l, TLorentzVector lead_p4, TLorentzVector sublead_p4, int mu_ind, float bdtoutput, float evweight, std::string label){

    if(mu_ind<0) {
        std::cout<<"mu_ind is "<<mu_ind<<std::endl;
        std::cout<<"leaving ControlPlotsMuonTag2012B"<<std::endl;
        return;
    }

    TLorentzVector* mymu=(TLorentzVector*) l.mu_glo_p4->At(mu_ind);
    int muVtx=l.FindMuonVertex(mu_ind);
    int mu_cat=(int)(abs(lead_p4.Eta())>1.5 || abs(sublead_p4.Eta())>1.5);

    l.FillHist(Form("MuonTag_mupt_%s",label.c_str()),       mu_cat, mymu->Pt(), evweight);
    l.FillHist(Form("MuonTag_mueta_%s",label.c_str()),      mu_cat, mymu->Eta(), evweight);
    l.FillHist(Form("MuonTag_muphi_%s",label.c_str()),      mu_cat, mymu->Phi(), evweight);
    float thisiso=((l.mu_glo_nehadiso04[mu_ind]+l.mu_glo_photiso04[mu_ind])>l.mu_dbCorr[mu_ind]) ?
    l.mu_glo_chhadiso04[mu_ind]+l.mu_glo_nehadiso04[mu_ind]+l.mu_glo_photiso04[mu_ind]-l.mu_dbCorr[mu_ind] : l.mu_glo_chhadiso04[mu_ind];
    l.FillHist(Form("MuonTag_muisoopt_%s",label.c_str()),   mu_cat, thisiso/mymu->Pt(), evweight);
    l.FillHist(Form("MuonTag_drmulead_%s",label.c_str()),   mu_cat, mymu->DeltaR(lead_p4), evweight);
    l.FillHist(Form("MuonTag_drmusub_%s",label.c_str()),    mu_cat, mymu->DeltaR(sublead_p4), evweight);
    l.FillHist(Form("MuonTag_d0_%s",label.c_str()),         mu_cat, l.mu_glo_D0Vtx[mu_ind][muVtx], evweight);
    l.FillHist(Form("MuonTag_dZ_%s",label.c_str()),         mu_cat, l.mu_glo_DZVtx[mu_ind][muVtx], evweight);
    l.FillHist(Form("MuonTag_diphomva_%s",label.c_str()),   mu_cat, bdtoutput, evweight);


}


bool PhotonAnalysis::VBFTag2011(LoopAll& l, int diphoton_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphoton_id==-1) return tag;

    TLorentzVector lead_p4    = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);

    std::pair<int, int> jets = l.Select2HighestPtJets(lead_p4, sublead_p4 );
    if(jets.first==-1 or jets.second==-1) return tag;

    TLorentzVector diphoton = lead_p4+sublead_p4;

    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.first);
    TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.second);
    TLorentzVector dijet = (*jet1) + (*jet2);

    myVBFLeadJPt= jet1->Pt();
    myVBFSubJPt = jet2->Pt();
    myVBF_Mjj   = dijet.M();
    myVBFdEta   = fabs(jet1->Eta() - jet2->Eta());
    myVBFZep    = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
    myVBFdPhi   = fabs(diphoton.DeltaPhi(dijet));
    myVBF_Mgg   = diphoton.M();

    if(nm1){
        tag = l.ApplyCutsFill(0, 1, eventweight, myweight);
    } else {
        tag = l.ApplyCuts(0, 1);
    }

    return tag;
}

bool PhotonAnalysis::VBFTag2012(int & ijet1, int & ijet2,
                LoopAll& l, int diphoton_id, float* smeared_pho_energy, bool nm1, float eventweight,
                float myweight, bool * jetid_flags)
{
    bool tag = false;
    bool getAngles = true;
    
    bool jetsPreselected=FillDijetVariables(ijet1, ijet2, l, diphoton_id, &smeared_pho_energy[0], jetid_flags, getAngles);
    
    if(PADEBUG) std::cout<<"VBFTag2012 jetsPreselected "<<jetsPreselected<<std::endl;
    if(jetsPreselected==false) return tag;
    
    if( mvaVbfSelection ) {
        if( myVBFLeadJPt>30. && myVBFSubJPt>20. && myVBF_Mjj > 250. ) { // FIXME hardcoded pre-selection thresholds
            if(nm1 && myVBF_Mgg>massMin && myVBF_Mgg<massMax) {
                l.FillCutPlots(0,1,"_nminus1",eventweight,myweight);
            }
            if (!multiclassVbfSelection || vbfVsDiphoVbfSelection ){
                myVBF_MVA = tmvaVbfReader_->EvaluateMVA(mvaVbfMethod);
                tag       = (myVBF_MVA > mvaVbfCatBoundaries.back());
            }
            else {
                myVBF_MVA0 = tmvaVbfReader_->EvaluateMulticlass(mvaVbfMethod)[0]; // signal vbf
                myVBF_MVA1 = tmvaVbfReader_->EvaluateMulticlass(mvaVbfMethod)[1]; // dipho
                myVBF_MVA2 = tmvaVbfReader_->EvaluateMulticlass(mvaVbfMethod)[2]; // gluglu
                //do transformation for bkg mvas
                myVBF_MVA1 = -myVBF_MVA1+1;
                myVBF_MVA2 = -myVBF_MVA2+1;
                tag        = (myVBF_MVA0 > multiclassVbfCatBoundaries0.back() && myVBF_MVA1 > multiclassVbfCatBoundaries1.back() && myVBF_MVA2 > multiclassVbfCatBoundaries2.back());
            }
            
            // this is moved to StatAnalysis::fillControlPlots
            // 	    if(nm1 && tag && myVBF_Mgg>massMin && myVBF_Mgg<massMax ) {
            // 		l.FillCutPlots(0,1,"_sequential",eventweight,myweight);
            // 	    }
            if( doDiphoMvaUpFront ) {
              if( vbfVsDiphoVbfSelection ) {
                tag = tag && ( l.dipho_BDT[diphoton_id] > multiclassVbfCatBoundaries1.back() );
              } else {
                tag = tag && ( l.dipho_BDT[diphoton_id] > bdtCategoryBoundaries.back() );
              }
            }
        }
    } else {
        if(nm1){
            tag = l.ApplyCutsFill(0,1, eventweight, myweight);
        } else {
            tag = l.ApplyCuts(0,1);
        }
    }
    
    if( mvaVbfSpin && (mvaVbfSelection || multiclassVbfSelection) )
    {
        myVBFSpin_Discriminant = tmvaVbfSpinReader_->EvaluateMVA(mvaVbfSpinMethod);
    }
    
    return tag;
}

bool PhotonAnalysis::VBFTag2013(int & ijet1, int & ijet2, LoopAll& l, int& diphotonVBF_id, float* smeared_pho_energy, 
                                bool vetodipho, bool kinonly, bool mvaselection, float eventweight, float myweight){
    bool tag = false;
    bool getAngles = false;
    
    diphotonVBF_id = l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVBFCut,subleadEtVBFCut,phoidMvaCut,applyPtoverM, 
                                                &smeared_pho_energy[0], vetodipho, kinonly );
    if(PADEBUG) std::cout<<"VBFTag2013 diphotonVBF_id "<<diphotonVBF_id<<std::endl;
    if(diphotonVBF_id==-1){
        return tag;
    } 
    
    bool jetsPreselected=FillDijetVariables(ijet1, ijet2, l, diphotonVBF_id, &smeared_pho_energy[0]);
    
    if(PADEBUG) std::cout<<"VBFTag2013 jetsPreselected "<<jetsPreselected<<std::endl;
    if(jetsPreselected==false) return tag;
    

    if(mvaselection){
        if( useGbrVbfMva ) {
            assert( doDiphoMvaUpFront );
            // Load di-photon MVA inputs for this di-photon so that we can calculate the 
            // combined BDT
            l.tmva_dipho_MIT_buf = l.tmva_dipho_MIT_cache.find(diphotonVBF_id)->second;
        }
        
        if(PADEBUG) std::cout<<"myVBFLeadJPt myVBFSubJPt myVBF_Mjj "<<myVBFLeadJPt<<" "<<myVBFSubJPt<<" "<<myVBF_Mjj<<std::endl;
        if( !(myVBFLeadJPt>30. && myVBFSubJPt>20. && myVBF_Mjj > 250.) ) { // FIXME hardcoded pre-selection thresholds
            return tag;
        }

        if( myVBF_Mgg>massMin && myVBF_Mgg<massMax) {
            l.FillCutPlots(0,1,"_nminus1",eventweight,myweight);
        }
        
        int vbfcat=-1;
        myVBFDIPHObdt   = l.dipho_BDT[diphotonVBF_id];
        myVBF_MVA       = (useGbrVbfMva ? gbrVbfReader_->eval()      : tmvaVbfReader_->EvaluateMVA(mvaVbfMethod)           );
        myVBFcombined   = (useGbrVbfMva ? gbrVbfDiphoReader_->eval() : tmvaVbfDiphoReader_->EvaluateMVA(mvaVbfDiphoMethod) );

        if(PADEBUG) std::cout<<"dipho dijet pt/m combined "<<myVBFDIPHObdt<<" "<<myVBF_MVA<<" "<<myVBFDiPhoPtOverM<<" "<<myVBFcombined<<std::endl;

        vbfcat=categoryFromBoundaries2D(multiclassVbfCatBoundaries0,multiclassVbfCatBoundaries1,multiclassVbfCatBoundaries2,
                                        myVBF_MVA,                  myVBFcombined,              -2);
       
        ///// std::cout << "VBFTag2013 " << l.run << " " << l.event << " " << myVBFLeadJEta << " " <<  myVBFSubJEta << " " << myVBFLeadJPt << " " << myVBFSubJPt << " " << myVBFZep << " " << myVBFdPhiTrunc << " " << myVBF_Mjj << " " << myVBFDiPhoPtOverM << " " << myVBFDIPHObdt << " " << myVBF_MVA << " " << myVBFcombined<<std::endl; 
 
        if(PADEBUG) std::cout<<"vbfcat dijet combinedmva "<<vbfcat<<" "<<myVBF_MVA<<" "<<myVBFcombined<<std::endl;
        if( vbfcat!=-1 ) tag = true;
    }
    return tag;
}


int PhotonAnalysis::categoryFromBoundaries(std::vector<float> & v, float val)
{
    int cat=-1;
    // val == v[0] would be -1; needs special condition
    if( val == v[0] ) { cat=0; }
    else {
        // bound is pointer to the ith boundary in v such that val>v[i]
        std::vector<float>::iterator bound =  lower_bound( v.begin(), v.end(), val, std::greater<float>  ());
        cat = ( val >= *bound ? bound - v.begin() - 1 : bound - v.begin() );
        if( cat >= v.size() - 1 ) { cat = -1; }
    }
    return cat;
}

int PhotonAnalysis::categoryFromBoundaries2D(std::vector<float> & v1, std::vector<float> & v2, std::vector<float> & v3, float val1, float val2, float val3 )
{
    int cat1temp =  categoryFromBoundaries(v1,val1);
    int cat2temp =  categoryFromBoundaries(v2,val2);
    int cat3temp =  categoryFromBoundaries(v3,val3);
    std::vector<int> vcat;
    vcat.push_back(cat1temp);
    vcat.push_back(cat2temp);
    vcat.push_back(cat3temp);
    int cat =  *min_element(vcat.begin(), vcat.end())==-1 ? -1 : *max_element(vcat.begin(), vcat.end());
    return cat;
}
bool PhotonAnalysis::FillDijetVariables(int & ijet1, int & ijet2, LoopAll& l, int diphoton_id, float* smeared_pho_energy,
                                        bool * jetid_flags, bool getAngles){
    bool filled=false;

    if(diphoton_id==-1) return filled;
    
    static std::vector<unsigned char> id_flags;

    if( jetid_flags == 0 ) {
        if(PADEBUG) std::cout<<"FillDijetVariable -- no id flags, re-making"<<std::endl;
        switchJetIdVertex( l, l.dipho_vtxind[diphoton_id] );
        id_flags.resize(l.jet_algoPF1_n);
        for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
            id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
        }
        jetid_flags = (bool*)&id_flags[0];
    }

    TLorentzVector lead_p4    = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
    if(PADEBUG) std::cout<<"FillDijetVariable -- photon ind "<<l.dipho_leadind[diphoton_id]<<"\t"<<l.dipho_subleadind[diphoton_id]<<std::endl;
    if(PADEBUG) std::cout<<"FillDijetVariable -- photon pt  "<<lead_p4.Pt()<<" "<<sublead_p4.Pt()<<std::endl;

    std::pair<int, int> jets;
    if(PADEBUG) std::cout<<"FillDijetVariable -- getting highest pt jets -- with PU jetveto?"<<usePUjetveto<<std::endl;
    if(usePUjetveto){
        jets = l.Select2HighestPtJets(lead_p4, sublead_p4, jetid_flags );
    } else {
        jets = l.Select2HighestPtJets(lead_p4, sublead_p4);
    }

    if(jets.first==-1 || jets.second==-1) {
        if(PADEBUG) std::cout<<"FillDijetVariable -- no jets"<<std::endl;
        return filled;
    }

    TLorentzVector diphoton = lead_p4+sublead_p4;

    ijet1 = jets.first; ijet2 = jets.second;
    if(PADEBUG) std::cout<<"FillDijetVariable -- ijet1 ijet2 "<<ijet1<<" "<<ijet2<<std::endl;
    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.first);
    TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.second);
    TLorentzVector dijet = (*jet1) + (*jet2);
    if(jet1->Pt() < jet2->Pt())
      std::swap(jet1, jet2);


    myVBFLeadJPt= jet1->Pt();
    myVBFSubJPt = jet2->Pt();
    myVBF_Mjj   = dijet.M();
    myVBFdEta   = fabs(jet1->Eta() - jet2->Eta());
    myVBFZep    = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
    myVBFdPhi   = fabs(diphoton.DeltaPhi(dijet));
    //    myVBFdPhiTrunc   = TMath::Min( (double)myVBFdPhi, TMath::Pi() - 0.2 );
    if(l.sqrtS==7){
        myVBFdPhiTrunc   = TMath::Min( (double)myVBFdPhi, 2.9416);
    } else if(l.sqrtS==8){
        myVBFdPhiTrunc   = TMath::Min( (double)myVBFdPhi, 2.916);
    }
    myVBF_Mgg   = diphoton.M();
    myVBFDiPhoPtOverM   = diphoton.Pt()   / myVBF_Mgg;
    myVBFLeadPhoPtOverM = lead_p4.Pt()    / myVBF_Mgg;
    myVBFSubPhoPtOverM  = sublead_p4.Pt() / myVBF_Mgg;
    myVBF_MVA  = -2.;
    myVBF_MVA0 = -2.;
    myVBF_MVA1 = -2.;
    myVBF_MVA2 = -2.;
    myVBFcombined = -2.;
    myVBFSpin_Discriminant = -2.;
    myVBF_deltaPhiGamGam = lead_p4.DeltaPhi(sublead_p4);
    myVBF_etaJJ = (jet1->Eta() + jet2->Eta())/2;
    myVBFLeadJEta = jet1->Eta();
    myVBFSubJEta = jet2->Eta();
    if(getAngles) VBFAngles(lead_p4, sublead_p4, *jet1, *jet2);

    filled=true;
    return filled;
}

bool PhotonAnalysis::VHhadronicTag2011(LoopAll& l, int& diphotonVHhad_id, float* smeared_pho_energy, bool *jetid_flags, bool mvaselection,bool vetodipho,bool kinonly){
    //francesco 

    bool tag = false;

    if(PADEBUG)    std::cout<<"-------zeds dead baby"<<std::endl;


    if(!mvaselection){
        diphotonVHhad_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHhadCut,subleadEtVHhadCut, 4,
                                                   applyPtoverM, &smeared_pho_energy[0], true);
    }else{
        diphotonVHhad_id=l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVHhadCut,subleadEtVHhadCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0],vetodipho,kinonly,
                                                   diphobdt_output_Cut_VHhad);
    }

    if(diphotonVHhad_id==-1) return tag;

    static std::vector<unsigned char> id_flags;
    if( jetid_flags == 0 ) {
        switchJetIdVertex( l, l.dipho_vtxind[diphotonVHhad_id] );
        id_flags.resize(l.jet_algoPF1_n);
        for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
            id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
        }
        jetid_flags = (bool*)&id_flags[0];
    }
    


    //////////////////Defining VH selection///////////////
    float ptLead_thresh,ptSublead_thresh,ptDiphot_thresh,ptLeadTrig_thresh,ptSubleadTrig_thresh;
    int nJets_thresh,nJetsUpper_thresh;
    float ptJets_thresh,mjjLower_thresh,mjjUpper_thresh,absCosThetaStar_thresh;

    //defining VH variables
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHhad_id], l.dipho_vtxind[diphotonVHhad_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHhad_id], l.dipho_vtxind[diphotonVHhad_id], &smeared_pho_energy[0]);

    std::pair<int, int> jets;
    if(usePUjetveto){
        jets = l.Select2HighestPtJets(lead_p4, sublead_p4, jetid_flags );
    } else {
        jets = l.Select2HighestPtJets(lead_p4, sublead_p4);
    }

    if(jets.first==-1 || jets.second==-1) return tag;



    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.first);
    TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.second);

    TLorentzVector dijet = (*jet1) + (*jet2);
    TLorentzVector diphoton = lead_p4+sublead_p4;

    int njets=0;
    int njets_looseptcut=0;
    int njets_btagloose=0;
    int njets_btagmedium=0;


    //Vh 0 btag selection
    //photon cuts
    ptLead_thresh=60.*diphoton.M()/120.;
    ptSublead_thresh=25.*diphoton.M()/120.;
    ptLeadTrig_thresh=33.;
    ptSubleadTrig_thresh=25.;

    ptDiphot_thresh=ptgg_0tag_cut*diphoton.M()/120.;//0 tag values
    absCosThetaStar_thresh=costhetastar_0tag_cut;


    //    absCosThetaStar_thresh=0.57;

    //jet cuts
    nJets_thresh=2;
    nJetsUpper_thresh=3;//no upper cut on 0 tag category
    ptJets_thresh=ptjet_0tag_cut;
    mjjLower_thresh=60.;
    mjjUpper_thresh=120.;

    //jet selection
    for(int ii=0; ii<l.jet_algoPF1_n; ++ii) {
	TLorentzVector * p4_jet = (TLorentzVector *) l.jet_algoPF1_p4->At(ii);
        if(jetid_flags != 0 && !jetid_flags[ii]) continue; 
        if(fabs(p4_jet->Eta()) > 2.4) continue;

        bool isJet_LeadPho = false;
	bool isJet_SubLeadPho = false;

	double dR_jet_PhoLead = p4_jet->DeltaR(lead_p4);
        if( dR_jet_PhoLead<0.5 ) isJet_LeadPho = true;

        double dR_jet_PhoSubLead = p4_jet->DeltaR(sublead_p4);
        if( dR_jet_PhoSubLead<0.5 ) isJet_SubLeadPho = true;

        if( isJet_LeadPho || isJet_SubLeadPho ) continue;

	if(p4_jet->Pt()>20.){
	    if(l.jet_algoPF1_csvBtag[ii]>0.244)njets_btagloose++;
	    if(l.jet_algoPF1_csvBtag[ii]>0.679)njets_btagmedium++;
	}

	if(p4_jet->Pt()>ptjet_loosecut)	njets_looseptcut++;

	if(p4_jet->Pt()<ptJets_thresh) continue;




	njets++;


	if(PADEBUG)std::cout<<"pt: "<<p4_jet->Pt()<<" btag_loose "<<njets_btagloose<<" btag_medium "<<njets_btagmedium<<std::endl;

    }


    //costhetaStar 

    TLorentzVector Vstar = dijet + diphoton;

    TLorentzVector H_Vstar(diphoton);
    H_Vstar.Boost(-Vstar.BoostVector());

    float cosThetaStar = -H_Vstar.CosTheta();
    float abs_cosThetaStar = fabs(cosThetaStar);



    bool isNotBtaggedLoose=(njets_btagloose==0);
    //    bool isEBEB=TMath::Abs(lead_p4.Eta())<1.4442 && TMath::Abs(sublead_p4.Eta())<1.4442;

    if(PADEBUG) std::cout<<" njets: "<<njets<<" abs_cosThetaStar: "<<abs_cosThetaStar<<std::endl;

    //doing the selection
    bool hasPassedCosThetaSelection=(abs_cosThetaStar<absCosThetaStar_thresh);
    bool hasPassedJetSelection= (njets>=nJets_thresh && njets_looseptcut<=nJetsUpper_thresh && isNotBtaggedLoose && dijet.M()>mjjLower_thresh && dijet.M()<mjjUpper_thresh);
    bool hasPassedPhotonSelection= (lead_p4.Pt()>ptLeadTrig_thresh && sublead_p4.Pt()> ptSubleadTrig_thresh &&lead_p4.Pt()> ptLead_thresh  &&  diphoton.Pt()>ptDiphot_thresh);



    if(hasPassedJetSelection && hasPassedPhotonSelection && hasPassedCosThetaSelection)tag=true;

      if (PADEBUG &&tag==true)  cout<<"tagged VH had, event"<<l.event<<"run "<<l.run<<" lumi "<<l.lumis<<endl;
    //   if(tag==true)  cout<<"tagged VH had, event"<<l.event<<"run "<<l.run<<" lumi "<<l.lumis<<endl;


    myVHhadLeadJPt = jet1->Pt();
    myVHhadSubJPt = jet2->Pt();
    myVHhad_Mjj = dijet.M();
    myVHhaddEta = fabs(jet1->Eta() - jet2->Eta());
    myVHhadZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
    myVHhaddPhi = fabs(diphoton.DeltaPhi(dijet));
    myVHhad_Mgg =diphoton.M();


    /*    if(nm1){
        tag = l.ApplyCutsFill(0,2, eventweight, myweight);
    } else {
        tag = l.ApplyCuts(0,2);
	}*/

    return tag;
}

bool PhotonAnalysis::VHhadronicTag2012(LoopAll& l, int& diphotonVHhad_id, float* smeared_pho_energy, bool *jetid_flags, bool mvaselection,bool vetodipho,bool kinonly){
    //one single category without btag requirements

    bool tag = false;


    if(!mvaselection){
        diphotonVHhad_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHhadCut,subleadEtVHhadCut, 4,
                                                   false, &smeared_pho_energy[0], true);
    }else{
        diphotonVHhad_id=l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVHhadCut,subleadEtVHhadCut,phoidMvaCut,
                                                   applyPtoverM, &smeared_pho_energy[0],vetodipho,kinonly,diphobdt_output_Cut_VHhad);
    }

    if(diphotonVHhad_id==-1) return tag;

    static std::vector<unsigned char> id_flags;
    if( jetid_flags == 0 ) {
        switchJetIdVertex( l, l.dipho_vtxind[diphotonVHhad_id] );
        id_flags.resize(l.jet_algoPF1_n);
        for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
            id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
        }
        jetid_flags = (bool*)&id_flags[0];
    }
    


    //////////////////Defining VH selection///////////////
    float ptDiphot_thresh;
    int nJets_thresh;
    float ptJets_thresh,mjjLower_thresh,mjjUpper_thresh,absCosThetaStar_thresh;

    //defining VH variables
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHhad_id], l.dipho_vtxind[diphotonVHhad_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHhad_id], l.dipho_vtxind[diphotonVHhad_id], &smeared_pho_energy[0]);

    std::pair<int, int> jets;
    if(usePUjetveto){
        jets = l.Select2HighestPtJets(lead_p4, sublead_p4, jetid_flags );
    } else {
        jets = l.Select2HighestPtJets(lead_p4, sublead_p4);
    }

    if(jets.first==-1 || jets.second==-1) return tag;


    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.first);
    TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.second);

    TLorentzVector dijet = (*jet1) + (*jet2);
    TLorentzVector diphoton = lead_p4+sublead_p4;

    int njets=0;
    int njets_looseptcut=0;
    int njets_btagloose=0;
    int njets_btagmedium=0;


    //Vh had selection
    //photon cuts
    ptDiphot_thresh=ptgg_0tag_cut*diphoton.M()/120.;//0 tag values
    absCosThetaStar_thresh=costhetastar_0tag_cut;



    //jet cuts
    nJets_thresh=2;
    ptJets_thresh=ptjet_0tag_cut;
    mjjLower_thresh=60.;
    mjjUpper_thresh=120.;

    //jet selection
    for(int ii=0; ii<l.jet_algoPF1_n; ++ii) {
        TLorentzVector * p4_jet = (TLorentzVector *) l.jet_algoPF1_p4->At(ii);
        if(jetid_flags != 0 && !jetid_flags[ii]) continue; 
        if(fabs(p4_jet->Eta()) > 2.4) continue;
        
        bool isJet_LeadPho = false;
        bool isJet_SubLeadPho = false;
        
        double dR_jet_PhoLead = p4_jet->DeltaR(lead_p4);
        if( dR_jet_PhoLead<0.5 ) isJet_LeadPho = true;
        
        double dR_jet_PhoSubLead = p4_jet->DeltaR(sublead_p4);
        if( dR_jet_PhoSubLead<0.5 ) isJet_SubLeadPho = true;
        
        if( isJet_LeadPho || isJet_SubLeadPho ) continue;
        
        if(p4_jet->Pt()<ptJets_thresh) continue;
        njets++;
                
        if(PADEBUG)std::cout<<"pt: "<<p4_jet->Pt()<<" btag_loose "<<njets_btagloose<<" btag_medium "<<njets_btagmedium<<std::endl;
    }


    //costhetaStar 

    TLorentzVector Vstar = dijet + diphoton;

    TLorentzVector H_Vstar(diphoton);
    H_Vstar.Boost(-Vstar.BoostVector());

    float cosThetaStar = -H_Vstar.CosTheta();
    float abs_cosThetaStar = fabs(cosThetaStar);



    if(PADEBUG) std::cout<<" njets: "<<njets<<" abs_cosThetaStar: "<<abs_cosThetaStar<<std::endl;

    //doing the selection
    bool hasPassedCosThetaSelection=(abs_cosThetaStar<absCosThetaStar_thresh);
    bool hasPassedJetSelection= (njets>=nJets_thresh && dijet.M()>mjjLower_thresh && dijet.M()<mjjUpper_thresh);
    bool hasPassedPhotonSelection=  diphoton.Pt()>ptDiphot_thresh;



    if(hasPassedJetSelection && hasPassedPhotonSelection && hasPassedCosThetaSelection)tag=true;

    if (PADEBUG &&tag==true)  cout<<"tagged VH had, event"<<l.event<<"run "<<l.run<<" lumi "<<l.lumis<<endl;

    myVHhadLeadJPt = jet1->Pt();
    myVHhadSubJPt = jet2->Pt();
    myVHhad_Mjj = dijet.M();
    myVHhaddEta = fabs(jet1->Eta() - jet2->Eta());
    myVHhadZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
    myVHhaddPhi = fabs(diphoton.DeltaPhi(dijet));
    myVHhad_Mgg =diphoton.M();

    return tag;
}



bool PhotonAnalysis::VHhadronicBtag2012(LoopAll& l, int& diphotonVHhadBtag_id, float* smeared_pho_energy, bool *jetid_flags, bool mvaselection,bool vetodipho,bool kinonly){
    //francesco 
    bool tag = false;



    if(!mvaselection){
        diphotonVHhadBtag_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHhadBtagCut,subleadEtVHhadBtagCut, 4,
                                                       applyPtoverM, &smeared_pho_energy[0], true);
    }else{
        diphotonVHhadBtag_id=l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVHhadBtagCut,subleadEtVHhadBtagCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0],
                                                       vetodipho,kinonly,diphobdt_output_Cut_VHhadBtag);
    }

    if(diphotonVHhadBtag_id==-1) return tag;

    static std::vector<unsigned char> id_flags;
    if( jetid_flags == 0 ) {
        switchJetIdVertex( l, l.dipho_vtxind[diphotonVHhadBtag_id] );
        id_flags.resize(l.jet_algoPF1_n);
        for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
            id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
        }
        jetid_flags = (bool*)&id_flags[0];
    }


    //////////////////Defining VH selection///////////////
    float ptLead_thresh,ptSublead_thresh,ptLeadTrig_thresh,ptSubleadTrig_thresh,ptDiphot_thresh;
    int nJets_thresh,nJetsUpper_thresh;
    float ptJets_thresh,mjjLower_thresh,mjjUpper_thresh,absCosThetaStar_thresh;

    //defining VH variables
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHhadBtag_id], l.dipho_vtxind[diphotonVHhadBtag_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHhadBtag_id], l.dipho_vtxind[diphotonVHhadBtag_id], &smeared_pho_energy[0]);


    std::pair<int, int> jets;
    if(usePUjetveto){
        jets = SelectBtaggedAndHighestPtJets(l,diphotonVHhadBtag_id,lead_p4, sublead_p4, jetid_flags );
    } else {
        jets = SelectBtaggedAndHighestPtJets(l,diphotonVHhadBtag_id,lead_p4, sublead_p4);
    }

    if(jets.first==-1 or jets.second==-1) return tag;

    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.first);
    TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.second);

    TLorentzVector dijet = (*jet1) + (*jet2);
    TLorentzVector diphoton = lead_p4+sublead_p4;

    int njets=0;
    int njets_looseptcut=0;
    int njets_btagloose=0;
    int njets_btagmedium=0;


    //Vh 1 btag selection
    //photon cuts
    ptLead_thresh=60.*diphoton.M()/120.;
    ptSublead_thresh=25.*diphoton.M()/120.;
    //    ptDiphot_thresh=70.*diphoton.M()/120.;//optimized values
    //    absCosThetaStar_thresh=0.84;
    ptDiphot_thresh=ptgg_btag_cut*diphoton.M()/120.;//0 tag values
    absCosThetaStar_thresh=costhetastar_btag_cut;


    //    cout<<"ptgg_btag_cut"<<ptgg_btag_cut<<endl;

    ptLeadTrig_thresh=33.;
    ptSubleadTrig_thresh=25.;


    //jet cuts
    nJets_thresh=2;
    nJetsUpper_thresh=3;
    //    ptJets_thresh=20.;//optimized values
    ptJets_thresh=ptjet_btag_cut;//0 tag values
    mjjLower_thresh=60.;
    mjjUpper_thresh=120.;


    //jet selection
    for(int ii=0; ii<l.jet_algoPF1_n; ++ii) {

	TLorentzVector * p4_jet = (TLorentzVector *) l.jet_algoPF1_p4->At(ii);
        if(jetid_flags != 0 && !jetid_flags[ii]) continue; 
        if(fabs(p4_jet->Eta()) > 2.4) continue;

        bool isJet_LeadPho = false;
	bool isJet_SubLeadPho = false;

	double dR_jet_PhoLead = p4_jet->DeltaR(lead_p4);
        if( dR_jet_PhoLead<0.5 ) isJet_LeadPho = true;

        double dR_jet_PhoSubLead = p4_jet->DeltaR(sublead_p4);
        if( dR_jet_PhoSubLead<0.5 ) isJet_SubLeadPho = true;

        if( isJet_LeadPho || isJet_SubLeadPho ) continue;

	if(ii==jets.first && l.jet_algoPF1_csvBtag[ii]>0.244){
	    njets++;
	    njets_btagloose++;
	}

	if(p4_jet->Pt()>ptjet_loosecut)	njets_looseptcut++;

	if(p4_jet->Pt()<ptJets_thresh) continue;


	njets++;

	if(l.jet_algoPF1_csvBtag[ii]>0.244)njets_btagloose++;
	if(l.jet_algoPF1_csvBtag[ii]>0.679)njets_btagmedium++;



		if(PADEBUG)std::cout<<"pt: "<<p4_jet->Pt()<<" btag_loose "<<njets_btagloose<<" btag_medium "<<njets_btagmedium<<std::endl;

    }



    //costhetaStar 

    TLorentzVector Vstar = dijet + diphoton;

    TLorentzVector H_Vstar(diphoton);
    H_Vstar.Boost(-Vstar.BoostVector());

    float cosThetaStar = -H_Vstar.CosTheta();
    float abs_cosThetaStar = fabs(cosThetaStar);


    bool isBtaggedLoose=(njets_btagloose>0);
    //    bool isEBEB=TMath::Abs(lead_p4.Eta())<1.4442 && TMath::Abs(sublead_p4.Eta())<1.4442;

       if(PADEBUG) std::cout<<" njets: "<<njets<<" abs_cosThetaStar: "<<abs_cosThetaStar<<std::endl;

    //doing the selection
    bool hasPassedCosThetaSelection=(abs_cosThetaStar<absCosThetaStar_thresh);
    bool hasPassedJetSelection= (njets>=nJets_thresh && njets_looseptcut<=nJetsUpper_thresh && isBtaggedLoose && dijet.M()>mjjLower_thresh && dijet.M()<mjjUpper_thresh);
    bool hasPassedPhotonSelection= (lead_p4.Pt()>ptLeadTrig_thresh && sublead_p4.Pt()> ptSubleadTrig_thresh && lead_p4.Pt()> ptLead_thresh  &&  diphoton.Pt()>ptDiphot_thresh);



    if(hasPassedJetSelection && hasPassedPhotonSelection && hasPassedCosThetaSelection)tag=true;

    if (PADEBUG && tag==true) cout<<"tagged VHBtag"<<endl;
    if(PADEBUG && tag==true)  cout<<"tagged VH had btag, event"<<l.event<<"run "<<l.run<<" lumi "<<l.lumis<<endl;

    if(PADEBUG && tag==true){
	std::cout<<"------------------------DEBUGGING------------------------"<<endl;
	std::cout<<"abs_cosThetaStar="<<abs_cosThetaStar<<" njets= "<<njets<<" nbtag loose="<<njets_btagloose<<" nbtag medium="<<njets_btagmedium
		 <<std::endl<<" ptphot1/2="<<lead_p4.Pt()<<"/"<<sublead_p4.Pt()<<" ptgg="<<diphoton.Pt()<<std::endl
		 <<" ptjet1="<<(*jet1).Pt()<<" ptjet2="<<(*jet2).Pt()<< " mjj"<<dijet.M()<<endl;
    }



    return tag;
}


bool PhotonAnalysis::TTHhadronicTag2012(LoopAll& l, int& diphotonTTHhad_id, float* smeared_pho_energy, bool *jetid_flags, bool mvaselection,bool vetodipho,bool kinonly){
    //francesco 
    bool tag = false;

    if(!mvaselection){
        diphotonTTHhad_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtTTHhadCut,subleadEtTTHhadCut, 4,
                                                    false, &smeared_pho_energy[0], true);
    }else{
        diphotonTTHhad_id=l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtTTHhadCut,subleadEtTTHhadCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0],
                                                    vetodipho,kinonly,diphobdt_output_Cut_TTHhad);
    }

    if(diphotonTTHhad_id==-1) return tag;

    static std::vector<unsigned char> id_flags;
    if( jetid_flags == 0 ) {
        switchJetIdVertex( l, l.dipho_vtxind[diphotonTTHhad_id] );
        id_flags.resize(l.jet_algoPF1_n);
        for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
            id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
        }
        jetid_flags = (bool*)&id_flags[0];
    }

    //////////////////Defining TTH selection///////////////
    float ptLead_thresh,ptSublead_thresh,ptLeadTrig_thresh,ptSubleadTrig_thresh;
    int nJets_thresh;
    float ptJets_thresh;

    //defining TTH variables
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonTTHhad_id], l.dipho_vtxind[diphotonTTHhad_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonTTHhad_id], l.dipho_vtxind[diphotonTTHhad_id], &smeared_pho_energy[0]);
    TLorentzVector diphoton = lead_p4+sublead_p4;


    int njets=0;
    int njets_btagloose=0;
    int njets_btagmedium=0;




    //photon cuts
    ptLead_thresh=60.*diphoton.M()/120.;
    ptSublead_thresh=25.*diphoton.M()/120.;
    ptLeadTrig_thresh=33.;
    ptSubleadTrig_thresh=25.;

    //jet cuts
    nJets_thresh=njets_tthHad_thresh;
    //    nJets_thresh=4;
    ptJets_thresh=ptJets_ttH_thresh;

    //jet selection
    for(int ii=0; ii<l.jet_algoPF1_n; ++ii) {

	TLorentzVector * p4_jet = (TLorentzVector *) l.jet_algoPF1_p4->At(ii);
        if(jetid_flags != 0 && !jetid_flags[ii]) continue; 
        if(fabs(p4_jet->Eta()) > 2.4) continue;

        bool isJet_LeadPho = false;
	bool isJet_SubLeadPho = false;

	double dR_jet_PhoLead = p4_jet->DeltaR(lead_p4);
        if( dR_jet_PhoLead<0.5 ) isJet_LeadPho = true;

        double dR_jet_PhoSubLead = p4_jet->DeltaR(sublead_p4);
        if( dR_jet_PhoSubLead<0.5 ) isJet_SubLeadPho = true;

        if( isJet_LeadPho || isJet_SubLeadPho ) continue;

	if(p4_jet->Pt()<ptJets_thresh) continue;


	njets++;

	if(l.jet_algoPF1_csvBtag[ii]>0.244)njets_btagloose++;
	if(l.jet_algoPF1_csvBtag[ii]>0.679)njets_btagmedium++;

	if(PADEBUG)
std::cout<<"pt: "<<p4_jet->Pt()<<" btag_loose "<<njets_btagloose<<" btag_medium "<<njets_btagmedium<<std::endl;




    }


    bool isBtaggedMedium;
    if(!removeBtagtth){
	isBtaggedMedium=(njets_btagmedium>0);
    }else{
	isBtaggedMedium=true;
    }
       if(PADEBUG)
	std::cout<<" njets: "<<njets<<std::endl;

    //doing the selection
    bool hasPassedJetSelection= (njets>=nJets_thresh && isBtaggedMedium);
    bool hasPassedPhotonSelection= (lead_p4.Pt()>ptLeadTrig_thresh && sublead_p4.Pt()> ptSubleadTrig_thresh && lead_p4.Pt()> ptLead_thresh);// && sublead_p4.Pt()>ptSublead_thresh);



    if(hasPassedJetSelection && hasPassedPhotonSelection)tag=true;

      if (PADEBUG && tag==true) cout<<"tagged TTH had"<<endl;

    //      if(tag==true)  cout<<"tagged TTHhad , event"<<l.event<<"run "<<l.run<<" lumi "<<l.lumis<<endl;


    return tag;
}

bool PhotonAnalysis::TTHTag7TeV(LoopAll& l, int& diphotonTTHlep_id, float* smeared_pho_energy, bool *jetid_flags, bool mvaselection,bool vetodipho,bool kinonly){
    bool tag=false;
    bool leptag=false;
    bool hadtag=false;
    if(!mvaselection){
        leptag=TTHleptonicTag2012(l, diphotonTTHlep_id, &smeared_pho_energy[0]);
        if(!leptag) hadtag=TTHhadronicTag2012(l, diphotonTTHlep_id, &smeared_pho_energy[0]);
    }else{
        leptag=TTHleptonicTag2012(l, diphotonTTHlep_id, &smeared_pho_energy[0], 0, true, vetodipho, kinonly);
        if(!leptag) hadtag=TTHhadronicTag2012(l, diphotonTTHlep_id, &smeared_pho_energy[0], 0, true, vetodipho, kinonly); 
    }
    if(leptag || hadtag)tag=true;
    return tag;
}

bool PhotonAnalysis::TTHleptonicTag2012(LoopAll& l, int& diphotonTTHlep_id, float* smeared_pho_energy, bool *jetid_flags, bool mvaselection,bool vetodipho,bool kinonly){
    //francesco 
    bool tag = false;

    int isLep_mu=0;
    int isLep_ele=0;
    int el_ind=-1;
    int mu_ind=-1;

    if(PADEBUG)
	std::cout<<"----------------this is tth lep"<<std::endl;


    //lepton requirement
    //defining TTH variables

    float myptcut=20.;
    int elInd = l.ElectronSelectionMVA2012(myptcut);
    int muonInd = l.MuonSelection2012B(myptcut);


    TLorentzVector* el_tag;
    TLorentzVector* mu_tag;

    bool passElePhotonCuts=false;
    bool passMuPhotonCuts=false;

    if(elInd != -1){
        el_tag = (TLorentzVector*) l.el_std_p4->At(elInd);
    }

    int elVtx=-1;
    std::vector<bool> veto_indices;
    veto_indices.clear();

    if(elInd!=-1) {
	TLorentzVector* myel = (TLorentzVector*) l.el_std_p4->At(elInd);
	TLorentzVector* myelsc = (TLorentzVector*) l.el_std_sc->At(elInd);

    float drtoveto = drSC_lep;
    float drgsftoveto = drGsf_lep;

    l.PhotonsToVeto(myelsc, drtoveto,veto_indices, true, drgsftoveto);
	elVtx=l.FindElectronVertex(elInd);

	// need to check again for d0 and dZ (couldn't before because we didn't have the vertex)                                        
    if(!(l.ElectronMVACuts(elInd, elVtx)))elInd=-1;
    if(elInd>-1)passElePhotonCuts=true;
    }


    if(!mvaselection){
        diphotonTTHlep_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtTTHlepCut,subleadEtTTHlepCut, 4,
                                                    false, &smeared_pho_energy[0], true, -1, veto_indices);
    }else{
        diphotonTTHlep_id=l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtTTHlepCut,subleadEtTTHlepCut,phoidMvaCut,applyPtoverM, &smeared_pho_energy[0],
                                                    vetodipho,kinonly,diphobdt_output_Cut_TTHlep,-1,false, veto_indices );
    }

    if(diphotonTTHlep_id==-1) return tag;

    //defining TTH variables                                                                                                                                 
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonTTHlep_id], l.dipho_vtxind[diphotonTTHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonTTHlep_id], l.dipho_vtxind[diphotonTTHlep_id], &smeared_pho_energy[0]);
    TLorentzVector diphoton = lead_p4+sublead_p4;

    if(muonInd != -1 && diphotonTTHlep_id !=1){
	mu_tag= (TLorentzVector*) l.mu_glo_p4->At(muonInd);
	passMuPhotonCuts=l.MuonPhotonCuts2012B(lead_p4, sublead_p4, mu_tag,drSC_lep);
    }

    if((elInd==-1) && (muonInd==-1))return tag;
    if(passElePhotonCuts == false && passMuPhotonCuts == false)return tag;


    if(muonInd != -1 && elInd==-1){
	if(passMuPhotonCuts){
	    isLep_mu=1;
	    mu_ind=muonInd;
	    el_ind=-1;
	}
    }
    if(elInd !=- 1 && muonInd ==-1){
	if(passElePhotonCuts){
	    isLep_ele=1;
	    el_ind=elInd;
	    mu_ind=-1;
	}
    }


    if(muonInd != -1 && elInd != -1){
	if(passMuPhotonCuts && passElePhotonCuts){
	    if(el_tag->Pt()<mu_tag->Pt()){
		isLep_mu=1;
		mu_ind=muonInd;
		el_ind=-1;
	    }else{
		isLep_ele=1;
		el_ind=elInd;
		mu_ind=-1;
	    }
	}else if(passMuPhotonCuts && !passElePhotonCuts){
	    isLep_mu=1;
	    mu_ind=muonInd;
	    el_ind=-1;
	}else if(passElePhotonCuts && !passMuPhotonCuts){
	    isLep_ele=1;
	    el_ind=elInd;
	    mu_ind=-1;
	}
    }

    if(isLep_ele!=1 && isLep_mu !=1) return false;
   

    static std::vector<unsigned char> id_flags;
    if( jetid_flags == 0 ) {
	switchJetIdVertex( l, l.dipho_vtxind[diphotonTTHlep_id] );
	id_flags.resize(l.jet_algoPF1_n);
	for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
	    id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
	}

	jetid_flags = (bool*)&id_flags[0];
    }





    /*    if(doDrGsfTrackCut){
	int leadind=l.dipho_leadind[diphotonTTHlep_id];
	int subleadind=l.dipho_subleadind[diphotonTTHlep_id];    
	//	std::cout<<"dr:"<<l.pho_drtotk_25_99[leadind]<<" "<<l.pho_drtotk_25_99[leadind]<<std::endl;
	if(l.pho_drtotk_25_99[leadind]<1.0 || l.pho_drtotk_25_99[subleadind]<1.0)return tag;

	}*/
    //cut on gsf track for electrons
    /*    bool passedDeltaRcut=true;
    if(elInd !=-1){
	for(int ii=0; ii<l.gsf_tk_n; ++ii) {
	    TLorentzVector * gsf_tk = (TLorentzVector *) l.gsf_tk_p4->At(ii);
	    float deltaRGsf_lead=gsf_tk->DeltaR(lead_p4);
	    float deltaRGsf_subLead=gsf_tk->DeltaR(sublead_p4);
	    if(deltaRGsf_lead<1.0 || deltaRGsf_subLead <1.0)passedDeltaRcut=false;
	    
	}
    }

    if(passedDeltaRcut==false)return tag;
    */

        if(PADEBUG)
	std::cout<<"elInd:"<<elInd<<" muonInd:"<<muonInd<<endl;

    //////////////////Defining TTH selection///////////////
	float ptLead_thresh,ptSublead_thresh,ptLeadTrig_thresh,ptSubleadTrig_thresh;
	int nJets_thresh;
	float ptJets_thresh;


    int njets=0;
    int njets_btagloose=0;
    int njets_btagmedium=0;


    //photon cuts
    ptLead_thresh=60.*diphoton.M()/120.;
    ptSublead_thresh=25.*diphoton.M()/120.;
    ptLeadTrig_thresh=33.;
    ptSubleadTrig_thresh=25.;

    //jet cuts
    nJets_thresh=2;
    ptJets_thresh=ptJets_ttH_thresh;

    TLorentzVector* lep;
    if(isLep_ele){
	lep= (TLorentzVector*) l.el_std_p4->At(elInd);
    }else if(isLep_mu){
	lep= (TLorentzVector*)l.mu_glo_p4->At(muonInd);
    }



    //jet selection
    for(int ii=0; ii<l.jet_algoPF1_n; ++ii) {

	TLorentzVector * p4_jet = (TLorentzVector *) l.jet_algoPF1_p4->At(ii);
        if(jetid_flags != 0 && !jetid_flags[ii]) continue; 
        if(fabs(p4_jet->Eta()) > 2.4) continue;
        bool isJet_LeadPho = false;
	bool isJet_SubLeadPho = false;
	bool isJet_Lep=false;

	double dR_jet_PhoLead = p4_jet->DeltaR(lead_p4);
        if( dR_jet_PhoLead<0.5 ) isJet_LeadPho = true;

        double dR_jet_PhoSubLead = p4_jet->DeltaR(sublead_p4);
        if( dR_jet_PhoSubLead<0.5 ) isJet_SubLeadPho = true;

	double dr_jet_lep= p4_jet->DeltaR(*lep);
	//	cout<<"-----"<<dr_jet_lep<<" ";
	if(dr_jet_lep<0.5) isJet_Lep = true;

        if( isJet_LeadPho || isJet_SubLeadPho || isJet_Lep ) continue;

	if(p4_jet->Pt()<ptJets_thresh) continue;

	//	cout<<ptJets_thresh<<endl;
	njets++;

	if(l.jet_algoPF1_csvBtag[ii]>0.244)njets_btagloose++;
	if(l.jet_algoPF1_csvBtag[ii]>0.679)njets_btagmedium++;

	if(PADEBUG)
	std::cout<<"pt: "<<p4_jet->Pt()<<" btag_loose "<<njets_btagloose<<" btag_medium "<<njets_btagmedium<<std::endl;

    }

    bool isBtaggedMedium;
    isBtaggedMedium=(njets_btagmedium>0);
    bool isBtaggedLoose;
    isBtaggedLoose=(njets_btagloose>0);	


    if(PADEBUG)
	std::cout<<" njets: "<<njets<<std::endl;

    //doing the selection
        bool hasPassedJetSelection= (njets>=nJets_thresh && isBtaggedMedium);
    //        bool hasPassedJetSelection = (njets>=nJets_thresh && isBtaggedLoose);
    bool hasPassedPhotonSelection= (lead_p4.Pt()>ptLeadTrig_thresh && sublead_p4.Pt()> ptSubleadTrig_thresh && lead_p4.Pt()> ptLead_thresh);// && sublead_p4.Pt()>ptSublead_thresh);



    if(hasPassedJetSelection && hasPassedPhotonSelection)tag=true;

    if (PADEBUG && tag==true) cout<<"tagged TTHlep"<<endl;
    //    if(tag==true)  cout<<"tagged TTHhad , event"<<l.event<<"run "<<l.run<<" lumi "<<l.lumis<<endl;
    
    return tag;
}


void PhotonAnalysis::computeBtagEff(LoopAll &l){
    for (int gi=0;gi<l.gp_n;gi++){
		if(l.gp_status[gi]!=3)continue;
		if (abs(l.gp_pdgid[gi])>5) continue;
		//		std::cout<<"gp:"<<l.gp_status[gi]<<" pdg"<<l.gp_status[gi]<<endl;
		TLorentzVector* jet_mc = (TLorentzVector*)l.gp_p4->At(gi);
		//std::cout<<"pt "<<jet_mc->Pt()<<std::endl;
		for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) {
		    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
		    if(fabs(jet1->Eta()) > 2.4) continue;
		    if(fabs(jet1->Pt()) < 25.) continue;
		    double dr_jet=jet1->DeltaR(*jet_mc);
		    //  std::cout<<"eta"<<fabs(jet1->Eta())<<" pt:"<<fabs(jet1->Pt())<<"dr_jet"<<dr_jet<<std::endl;
		    if(dr_jet>0.3) continue;//mc matching

		    //		    std::cout<<l.gp_pdgid[gi]<<" "<<l.jet_algoPF1_csvBtag[ijet]<<endl;
		    if(abs(l.gp_pdgid[gi])==5){
			nBMC++;
			if(l.jet_algoPF1_csvBtag[ijet]>0.244)nBtagB++;
		    }else if(abs(l.gp_pdgid[gi])==4){
			nCMC++;
			if(l.jet_algoPF1_csvBtag[ijet]>0.244)nBtagC++;
		    }else{
			nLMC++;
			if(l.jet_algoPF1_csvBtag[ijet]>0.244)nBtagL++;
		    }

		    //		    std::cout<<nBMC<<" "<<nCMC<<" "<<nBMC<<endl;

		}//ijet
    }//gp

    return;
}


float PhotonAnalysis::ElectronSFReweight(LoopAll &l){

    float scale_factor=-10;
    float el_ind=0;

    TLorentzVector* thissc = (TLorentzVector*) l.el_std_sc->At(el_ind);
    float thiseta = fabs(thissc->Eta());

    // central value                                                                                                                                                        
    if ( fabs(thiseta)<1. )                           scale_factor = 0.984;                                                                                               
    if ( fabs(thiseta)<1.442 && fabs(thiseta)>=1.0 )  scale_factor = 0.988;                                                                                               
    if ( fabs(thiseta)<2.1 && fabs(thiseta)>=1.566 )  scale_factor = 1.016;                                                                                               
    if ( fabs(thiseta)>=2.1 )                         scale_factor = 0.999;                                                                                               
    // + 1 sigma                                                                                                                                                            
    //if ( fabs(thiseta)<1. )                           scale_factor = 0.986;                                                                                               
    //if ( fabs(thiseta)<1.442 && fabs(thiseta)>=1.0 )  scale_factor = 0.999;                                                                                               
    //if ( fabs(thiseta)<2.1 && fabs(thiseta)>=1.566 )  scale_factor = 1.024;                                                                                               
    //if ( fabs(thiseta)>=2.1 )                         scale_factor = 1.017;                                                                                               
    //                                                                                                                                                                      
    // to switch off                                                                                                                                                        
    // scale_factor = 1.;        

    return scale_factor;

}



float PhotonAnalysis::MuonSFReweight(LoopAll &l){

   float scale_factor=-10;
   int mu_ind=0;

    TLorentzVector* mymu = (TLorentzVector*) l.mu_glo_p4->At(mu_ind);
    float etaMu = mymu->Eta();
    // values for full dataset of moriond https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=233592
    if ( fabs(etaMu)<0.9 )                     scale_factor = 0.9938;                                                                                                     
    if ( fabs(etaMu)>=0.9 && fabs(etaMu)<1.2 ) scale_factor = 0.9915;                                                                                                     
    if ( fabs(etaMu)>=1.2 )                    scale_factor = 0.9992;                                                                                                     
    // + 1 sigma                                                                                                                                                            
    //if ( fabs(etaMu)<0.9 )                     scale_factor = 0.9983;                                                                                                     
    //if ( fabs(etaMu)>=0.9 && fabs(etaMu)<1.2 ) scale_factor = 0.9979;                                                                                                     
    //if ( fabs(etaMu)>=1.2 )                    scale_factor = 1.0056;                                                                                                     
    // no SF                                                                                                                                                                
    //scale_factor = 1.;                                                                                                                                                    

    return scale_factor;

}






float PhotonAnalysis::BtagReweight(LoopAll &l, bool shiftBtagEffUp_bc, bool shiftBtagEffDown_bc, bool shiftBtagEffUp_l, bool shiftBtagEffDown_l, int WP){
	    int nBjets=0,nCjets=0,nLjets=0; 
	    float eff_b,eff_c,eff_l;
	    float SFb,SFc,SFl;	   

	    //loose WP
	    if(WP==0){
		eff_b=0.80,eff_c=0.39,eff_l=0.13;
		SFb=1.008,SFc=1.008,SFl=1.09;
	    } else if(WP ==1) {
	    //medium WP
		eff_b=0.675,eff_c=0.08,eff_l=0.016;
		SFb=0.963,SFc=0.963,SFl=1.16;
	    }

	    if(WP == 0){
		if(shiftBtagEffUp_bc){
		    SFb=1.031,SFc=1.054,SFl=1.09;
		}else if(shiftBtagEffDown_bc){
		    SFb=0.985,SFc=0.962,SFl=1.09;
		}else if(shiftBtagEffUp_l){
		    SFb=1.008,SFc=1.008,SFl=1.11;
		}else if(shiftBtagEffDown_l){
		    SFb=1.008,SFc=1.008,SFl=1.07;
		}
	    }else if(WP ==1){
		if(shiftBtagEffUp_bc){
		    SFb=0.983,SFc=1.003,SFl=1.16;
		}else if(shiftBtagEffDown_bc){
		    SFb=0.943,SFc=0.923,SFl=1.16;
		}else if(shiftBtagEffUp_l){
		    SFb=0.963,SFc=0.963,SFl=1.20;
		}else if(shiftBtagEffDown_l){
		    SFb=0.963,SFc=0.963,SFl=1.12;
		}
	    }

	    float w=-1;
	    //	    	    std::cout<<"gp:"<<l.gp_n<<endl;//decomment on analysis input!
	    for (int gi=0;gi<l.gp_n;gi++){
		if(l.gp_status[gi]!=3)continue;
		if (abs(l.gp_pdgid[gi])>5 && abs(l.gp_pdgid[gi])!=21 ) continue;
		TLorentzVector* jet_mc = (TLorentzVector*)l.gp_p4->At(gi);
		if (jet_mc->Pt()<0.1 ) continue;
		//		std::cout<<"pt "<<jet_mc->Pt()<<std::endl;

		for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) {
		    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
		    if(fabs(jet1->Eta()) > 2.4) continue;
		    if(fabs(jet1->Pt()) < 25.) continue;
		    //		    cout<<"jet1 preso"<<endl;
		    double dr_jet=jet1->DeltaR(*jet_mc);
		    if(dr_jet>0.3) continue;//mc matching
		    if(abs(l.gp_pdgid[gi])==5){
			nBjets++;
		    }else if(abs(l.gp_pdgid[gi])==4){
			nCjets++;
		    }else{
			nLjets++;
		    }
		    //		    if(l.gp_pdgid[gi]==9 || l.gp_pdgid[gi]==21)cout<<"pdg"<<l.gp_pdgid[gi]<<endl;
		}
	    }

	    //	    cout<<"nnnn"<<nBjets<<" "<<nCjets<<" "<<nLjets<<endl;
	
	    BTagWeight b(1,6,eff_b,eff_c,eff_l);//b(nbtagmin,nbtagmax,beff,ceff,leff)


	    w= b.weight(nBjets,nCjets,nLjets,SFb,SFc,SFl,1);
	    /*
	      float wbcp= b.weight(ib,ic,il,1,1,1.1,1);
	      float wbcm= b.weight(ib,ic,il,0.8,0.8,1.1,1);
	      float wlp= b.weight(ib,ic,il,0.9,0.9,1.2,1);
	      float wlm= b.weight(ib,ic,il,0.9,0.9,1.0,1);
	      float err=sqrt((wbcp-wbcm)*(wbcp-wbcm)+(wlp-wlm)*(wlp-wlm));
	      std::cout << ib << "   " << ic << "  " << "  " << il << "  " << w << " +- " <<   std::endl; */
	   
	    //	    cout<<w<<endl;

	    if(w>0){
		return w;
	    }else{
		return 1;
	    }
}

float PhotonAnalysis::BtagReweight2013(LoopAll &l, bool shiftBtagEffUp_bc, bool shiftBtagEffDown_bc, bool shiftBtagEffUp_l, bool shiftBtagEffDown_l, int WP){//same procedure as BtagReweighting but with different mc variables to save space
	    int nBjets=0,nCjets=0,nLjets=0; 
	    float eff_b,eff_c,eff_l;
	    float SFb,SFc,SFl;	   

	    //loose WP
	    if(WP==0){
		eff_b=0.80,eff_c=0.39,eff_l=0.13;
		SFb=1.008,SFc=1.008,SFl=1.09;
	    } else if(WP ==1) {
	    //medium WP
		eff_b=0.675,eff_c=0.08,eff_l=0.016;
		SFb=0.963,SFc=0.963,SFl=1.16;
	    }

	    if(WP == 0){
		if(shiftBtagEffUp_bc){
		    SFb=1.031,SFc=1.054,SFl=1.09;
		}else if(shiftBtagEffDown_bc){
		    SFb=0.985,SFc=0.962,SFl=1.09;
		}else if(shiftBtagEffUp_l){
		    SFb=1.008,SFc=1.008,SFl=1.11;
		}else if(shiftBtagEffDown_l){
		    SFb=1.008,SFc=1.008,SFl=1.07;
		}
	    }else if(WP ==1){
		if(shiftBtagEffUp_bc){
		    SFb=0.983,SFc=1.003,SFl=1.16;
		}else if(shiftBtagEffDown_bc){
		    SFb=0.943,SFc=0.923,SFl=1.16;
		}else if(shiftBtagEffUp_l){
		    SFb=0.963,SFc=0.963,SFl=1.20;
		}else if(shiftBtagEffDown_l){
		    SFb=0.963,SFc=0.963,SFl=1.12;
		}
	    }

	    float w=-1;
	    //	    	    std::cout<<"jet_algoPF1_cgenMatched:"<<l.gp_n<<endl;//decomment on analysis input!
		for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet) {
		    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(ijet);
		    if(fabs(jet1->Eta()) > 2.4) continue;
		    if(fabs(jet1->Pt()) < 25.) continue;
            if(l.jet_algoPF1_bgenMatched[ijet]){
			nBjets++;
		    }else if(l.jet_algoPF1_cgenMatched[ijet]){
			nCjets++;
		    }else if(l.jet_algoPF1_lgenMatched[ijet]){
			nLjets++;
		    }
		}

	    //	    cout<<"nnnn"<<nBjets<<" "<<nCjets<<" "<<nLjets<<endl;
	
	    BTagWeight b(1,6,eff_b,eff_c,eff_l);//b(nbtagmin,nbtagmax,beff,ceff,leff)

	    w= b.weight(nBjets,nCjets,nLjets,SFb,SFc,SFl,1);

	    if(w>0){
		return w;
	    }else{
		return 1;
	    }
}


// ----------------------------------------------------------------------------------------------------



//met at analysis step
bool PhotonAnalysis::METTag2012(LoopAll& l, int& diphotonVHmet_id, float* smeared_pho_energy){
    bool tag = false;
    diphotonVHmet_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHmetCut, subleadEtVHmetCut, 4, false, &smeared_pho_energy[0], true);
    if(diphotonVHmet_id>-1) {
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHmet_id], l.dipho_vtxind[diphotonVHmet_id] , &smeared_pho_energy[0]);
        TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHmet_id], l.dipho_vtxind[diphotonVHmet_id] , &smeared_pho_energy[0]);
        TLorentzVector TwoPhoton_Vector = (lead_p4) + (sublead_p4);
        float m_gamgam = TwoPhoton_Vector.M();

        MetCorrections2012_Simple( l, lead_p4 , sublead_p4 );

        if (m_gamgam>100 && m_gamgam<180) {
          met_sync << " run: " << l.run
            << "\tevent: " << l.event
            << "\tleadPt: " << lead_p4.Pt()
            << "\tsubleadPt: " << sublead_p4.Pt()
            << "\tdiphomass: " << m_gamgam
            << "\traw_met: " << l.met_pfmet
            << "\traw_met_phi: " << l.met_phi_pfmet
            << "\tshifted_met: " << l.shiftMET_pt
            << "\tcorrected_met: " << l.shiftscaleMET_pt
            << "\tcorrected_met_phi: " << l.shiftscaleMET_phi
            << "\tjet_algoPF1_n: " << l.jet_algoPF1_n
            << endl;
        }
    }

    if(diphotonVHmet_id==-1) return tag;
    if( l.correctedpfMET > 70 ) tag = true;

    return tag;

}


bool PhotonAnalysis::METTag2012B(LoopAll& l, int& diphotonVHmet_id, int& met_cat, float* smeared_pho_energy, ofstream& met_sync, bool mvaselection, float phoidMvaCut, bool useUncorrMet){
    bool tag = false;

    int metVtx=0;  // use default
    if(mvaselection) {
        diphotonVHmet_id = l.DiphotonMITPreSelection(bdtTrainingType.c_str(),leadEtVHmetCut,subleadEtVHmetCut,phoidMvaCut,
            applyPtoverM, &smeared_pho_energy[0], true, true,diphobdt_output_Cut_VHMet);
    } else {
        diphotonVHmet_id = l.DiphotonCiCSelection( l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtVHmetCut,subleadEtVHmetCut, 4,
                                                   applyPtoverM, &smeared_pho_energy[0], true, -1);
    }

    if(diphotonVHmet_id!=-1){
	//	std::cout << "+++PFMET UNCORR " << l.met_pfmet << std::endl;
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHmet_id], metVtx, &smeared_pho_energy[0]);
        TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHmet_id], metVtx, &smeared_pho_energy[0]);
        TLorentzVector dipho_p4 = lead_p4+sublead_p4;

        float mass = dipho_p4.M();
        met_cat=(int)(abs(lead_p4.Eta())>1.5 || abs(sublead_p4.Eta())>1.5);

        tag = l.METAnalysis2012B(lead_p4, sublead_p4, useUncorrMet, true, moriond2013MetCorrection);
        if (tag) {
            if(met_cat!=0) tag=false;
        }

        if(!tag) diphotonVHmet_id=-1;

        // for synchronization
        if(tag){
            met_sync<<"run="<<l.run<<"\t";
            met_sync<<"lumis"<<l.lumis<<"\t";
            met_sync<<"event="<<(unsigned int) l.event<<"\t";
            met_sync<<"pt1="<<lead_p4.Pt()<<"\t";
            met_sync<<"eta1="<<lead_p4.Eta()<<"\t";
            met_sync<<"pt2="<<sublead_p4.Pt()<<"\t";
            met_sync<<"eta2="<<sublead_p4.Eta()<<"\t";
            met_sync<<"ptgg="<<dipho_p4.Pt()<<"\n";
            // met_sync<<"met="<<MET<<"\n";
        }
    }

    return tag;
}


void PhotonAnalysis::reVertex(LoopAll & l)
{
    l.vtx_std_ranked_list->clear();
    l.vtx_std_evt_mva->clear();
    std::vector<int> preselAll;

    for(int i=0; i<l.vtx_std_n ; i++) {
	preselAll.push_back(i);
    }
    vtxAna_.preselection( preselAll );

    for(int id=0; id<l.dipho_n; ++id ) {

	vtxAna_.setPairID(id);

	if( rematchConversions ) {
	    vtxAna_.setNConv(0);
	    PhotonInfo p1 = l.fillPhotonInfos(l.dipho_leadind[id], vtxAlgoParams.useAllConversions, 0); // WARNING using default photon energy: it's ok because we only re-do the conversion part
	    PhotonInfo p2 = l.fillPhotonInfos(l.dipho_subleadind[id], vtxAlgoParams.useAllConversions, 0); // WARNING using default photon energy: it's ok because we only re-do the conversion part

	    float zconv=0., szconv=0.;
	    vtxAna_.getZFromConvPair(zconv,  szconv, p1, p2);

	    for(int vid=0; vid<l.vtx_std_n; ++vid) {
		if( vtxAna_.nconv(vid) > 0 ) {
		    assert( szconv > 0. );
		    vtxAna_.setPullToConv( vid, fabs(  ((TVector3 *)l.vtx_std_xyz->At(vid))->Z() - zconv ) / szconv );
		} else {
		    vtxAna_.setPullToConv( vid, -1. );
		}
	    }
	}

	l.vtx_std_ranked_list->push_back( vtxAna_.rank(*tmvaPerVtxReader_,tmvaPerVtxMethod) );
	l.dipho_vtxind[id] = l.vtx_std_ranked_list->back()[0];
	if( tmvaPerEvtReader_ ) {
	    float vtxEvtMva;
        if (rescaleDZforVtxMVA) vtxEvtMva = vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back(),targetsigma/sourcesigma );
        else vtxEvtMva = vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back());
	    l.vtx_std_evt_mva->push_back(vtxEvtMva);
	}
    }
}

float PhotonAnalysis::PtReweight(double pt, int cur_type){
   
   if (cur_type==-37 || cur_type==-38 || cur_type==-39 || cur_type==-40) {
        return ptreweighHistSM->GetBinContent(ptreweighHistSM->FindBin(pt));      
   }
   else if (cur_type==-137) {
        return ptreweighHistGG->GetBinContent(ptreweighHistGG->FindBin(pt));      
   }
   else if (cur_type==-138) {
        return ptreweighHistQQ->GetBinContent(ptreweighHistQQ->FindBin(pt));      
   }
   else return 1.;
}

float PhotonAnalysis::BeamspotReweight(double vtxZ, double genZ) {
    if (genZ<(-100)) return 1.0;

    float diffVar = vtxZ-genZ;
    if (TMath::Abs(diffVar)<0.1) return 1.0;

    // gets to 4.8
    float newBSmean1  = 9.9391e-02;
    float newBSmean2  = 1.8902e-01;
    float newBSnorm1  = 5.3210e+00;
    float newBSnorm2  = 4.1813e+01;
    float newBSsigma1 = 9.7530e-01;
    float newBSsigma2 = 7.0811e+00;

    // gets to 5
    //float newBSmean1  = 9.549e-02;
    //float newBSmean2  = 2.334e-01;
    //float newBSnorm1  = 5.067e+00;
    //float newBSnorm2  = 4.159e+01;
    //float newBSsigma1 = 9.498e-01;
    //float newBSsigma2 = 7.289e+00;

    float oldBSmean1  = 7.2055e-02;
    float oldBSmean2  = 4.9986e-01;
    float oldBSnorm1  = 3.5411e+00;
    float oldBSnorm2  = 4.0258e+01;
    float oldBSsigma1 = 7.9678e-01;
    float oldBSsigma2 = 8.5356e+00;

    float newBSgaus1 = newBSnorm1*exp(-0.5*pow((diffVar-newBSmean1)/newBSsigma1,2));
    float newBSgaus2 = newBSnorm2*exp(-0.5*pow((diffVar-newBSmean2)/newBSsigma2,2));
    float oldBSgaus1 = oldBSnorm1*exp(-0.5*pow((diffVar-oldBSmean1)/oldBSsigma1,2));
    float oldBSgaus2 = oldBSnorm2*exp(-0.5*pow((diffVar-oldBSmean2)/oldBSsigma2,2));

    float reweight =  1.1235 * (newBSgaus1+newBSgaus2)/(oldBSgaus1+oldBSgaus2);
    //float reweight =  1.10332 * (newBSgaus1+newBSgaus2)/(oldBSgaus1+oldBSgaus2);

    //float sourceweight = exp(-pow(hardInterZ-beamspotZ,2)/2.0/sourcesigma/sourcesigma)/sourcesigma;
    //float targetweight = exp(-pow(hardInterZ-beamspotZ,2)/2.0/targetsigma/targetsigma)/targetsigma;

    //float reweight = targetweight/sourceweight;

    //if(PADEBUG) std::cout<<"hardInterZ targetweight/sourceweight reweight "<<hardInterZ<<" "<<targetweight<<"/"<<sourceweight<<" "<<reweight<<std::endl;
    if(PADEBUG) std::cout<< "vtxZ genZ newBSgaus1 newBSgaus2 oldBSgaus1 oldBSgaus2 reweight "<< vtxZ << " " << genZ << " " << newBSgaus1 << " " << newBSgaus2 << " " << oldBSgaus1 << " " << oldBSgaus2 << " " << reweight << " " << std::endl;

    return reweight;
}

void PhotonAnalysis::saveBSTrees(LoopAll &l, float evweight, int category, TLorentzVector Higgs, TVector3 *chosenVtx, TVector3 *genVtx, float diphobdt_output){

  l.FillTree("weight",evweight,"bstrees");
  l.FillTree("category",category,"bstrees");
  l.FillTree("bdtoutput",diphobdt_output,"bstrees");
  l.FillTree("mass",Higgs.M(),"bstrees");
  l.FillTree("ZfromGenToChosen",(*chosenVtx-*genVtx).Z(),"bstrees");
  double distance=1000.;
  int index;
  for (int i=0; i<l.vtx_std_n; i++){
    TVector3 *tempVtx = (TVector3*)l.vtx_std_xyz->At(i);
    if (TMath::Abs(distance)>TMath::Abs((*tempVtx-*genVtx).Z())){
        distance = (*tempVtx-*genVtx).Z();
        index=i;
    }
  }
  l.FillTree("ZfromGenToBest",distance,"bstrees");

}

float PhotonAnalysis::ComputeEventScaleError(LoopAll& l, int ipho1, int ipho2, float & scale1, float & scale1_err, float & scale2, float & scale2_err) {
    PhotonReducedInfo pho1 = photonInfoCollection[ipho1];
    PhotonReducedInfo pho2 = photonInfoCollection[ipho2];
    
    std::string cat1=eScaleSmearer->photonCategory(pho1);
    std::string cat2=eScaleSmearer->photonCategory(pho2);
    if( PADEBUG ) std::cout<<"cat1 "<<cat1<<std::endl;
    if( PADEBUG ) std::cout<<"cat2 "<<cat2<<std::endl;

    scale1 = 1.0;
    scale2 = 1.0;

    scale1_err = eScaleSmearer->myParameters_.scale_offset_error[cat1];
    scale2_err = eScaleSmearer->myParameters_.scale_offset_error[cat2];
    if( PADEBUG ) std::cout<<"scale_syst_err1 "<<scale1_err<<std::endl;
    if( PADEBUG ) std::cout<<"scale_syst_err2 "<<scale2_err<<std::endl;

    float scale_event = -1;
    float scale_event_err = -1;

    if(cat1==cat2){
        scale_event=scale1;
        scale_event_err=scale1_err;
    } else {
        scale_event=sqrt(scale1*scale2);
        scale_event_err=0.5*scale_event*sqrt(
                (scale1_err*scale1_err)/(scale1*scale2) +
                (scale2_err*scale2_err)/(scale2*scale2) );
    }

    return scale_event_err;
}

float PhotonAnalysis::ComputeEventSmearError(LoopAll& l, int ipho1, int ipho2, float & smear1, float & smear1_err, float & smear2, float & smear2_err) {
    PhotonReducedInfo pho1 = photonInfoCollection[ipho1];
    PhotonReducedInfo pho2 = photonInfoCollection[ipho2];
    
    std::string cat1=eResolSmearer->photonCategory(pho1);
    std::string cat2=eResolSmearer->photonCategory(pho2);
    if( PADEBUG ) std::cout<<"cat1 "<<cat1<<std::endl;
    if( PADEBUG ) std::cout<<"cat2 "<<cat2<<std::endl;

    smear1 = EnergySmearer::getSmearingSigma(eResolSmearer->myParameters_,cat1,pho1.energy(),pho1.caloPosition().Eta(),0.);
    smear2 = EnergySmearer::getSmearingSigma(eResolSmearer->myParameters_,cat2,pho2.energy(),pho2.caloPosition().Eta(),0.);
    if( PADEBUG ) std::cout<<"smear1 "<<smear1<<std::endl;
    if( PADEBUG ) std::cout<<"smear2 "<<smear2<<std::endl;

    smear1_err = eResolSmearer->myParameters_.smearing_sigma_error[cat1];
    smear2_err = eResolSmearer->myParameters_.smearing_sigma_error[cat2];
    if( PADEBUG ) std::cout<<"smear_stat_err1 "<<smear1_err<<std::endl;
    if( PADEBUG ) std::cout<<"smear_stat_err2 "<<smear2_err<<std::endl;

    float smear_event = -1;
    float smear_event_err = -1;

    if(cat1==cat2) {
        smear_event=sqrt(0.5)*smear1;
        smear_event_err=sqrt(0.5)*smear1_err;
    } else {
        smear_event=0.5*sqrt(smear1*smear1 + smear2*smear2);
        smear_event_err=0.25/smear_event*sqrt(
                (smear1*smear1*smear1_err*smear1_err) +
                (smear2*smear2*smear2_err*smear2_err) );
    }

    return smear_event_err;
}

pair<double,double> PhotonAnalysis::ComputeNewSigmaMs(LoopAll &l, int ipho1, int ipho2, int ivtx, float sys_shift){

    PhotonReducedInfo pho1 = photonInfoCollection[ipho1];
    PhotonReducedInfo pho2 = photonInfoCollection[ipho2];

    pho1.setCorrEnergyErr(pho1.corrEnergyErr()*(1.+sys_shift*0.1));
    pho2.setCorrEnergyErr(pho2.corrEnergyErr()*(1.+sys_shift*0.1));

    MassResolution tempMassRes;

    tempMassRes.Setup(l,&pho1,&pho2,ivtx,massResoPars, nR9Categories, nEtaCategories,beamspotSigma,true);
    double sigMright = tempMassRes.massResolutionEonlyNoSmear();
    double sigMwrong = tempMassRes.massResolutionWrongVtxNoSmear();
    pair<double,double> result(sigMright,sigMwrong);
    return result;
}


void PhotonAnalysis::saveDatCardTree(LoopAll &l, int cur_type, int category, int inc_cat, float evweight, int ipho1, int ipho2, int ivtx, TLorentzVector lead_p4, TLorentzVector sublead_p4, bool isCutBased, string proc, double sigmaMrv, double sigmaMwv, double sigmaMeonly, float vtxP, string trainPhil, string bdtType, float lead_id_mva, float sublead_id_mva){

   // track the scale and smear uncertainties per event
    float scale1, scale1_err, scale2, scale2_err;
    float smear1, smear1_err, smear2, smear2_err;
    float scale_err = ComputeEventScaleError(l,ipho1,ipho2,scale1,scale1_err,scale1,scale1_err);
    float smear_err = ComputeEventSmearError(l,ipho1,ipho2,smear1,smear1_err,smear1,smear1_err);

    if (!isCutBased){ 
        float bdtout = l.diphotonMVA(-1,ipho1,ipho2,ivtx,vtxP,lead_p4,sublead_p4,sigmaMrv, sigmaMwv, sigmaMeonly, trainPhil.c_str(), bdtType.c_str(), lead_id_mva, sublead_id_mva);
        
        // calculate diphobdt given shift in idMVA
        float bdtout_id_up   = l.diphotonMVA(-1,ipho1,ipho2,ivtx,vtxP,lead_p4,sublead_p4,sigmaMrv, sigmaMwv, sigmaMeonly, trainPhil.c_str(), bdtType.c_str(), lead_id_mva+0.01, sublead_id_mva+0.01);
        float bdtout_id_down = l.diphotonMVA(-1,ipho1,ipho2,ivtx,vtxP,lead_p4,sublead_p4,sigmaMrv, sigmaMwv, sigmaMeonly, trainPhil.c_str(), bdtType.c_str(), lead_id_mva-0.01, sublead_id_mva-0.01);
        
        // calculate diphobdt given shift in sigmaE from regression
        pair<double,double> newSigmaMsUp = ComputeNewSigmaMs(l,ipho1,ipho2,ivtx,1.);
        pair<double,double> newSigmaMsDown = ComputeNewSigmaMs(l,ipho1,ipho2,ivtx,-1.);
        float bdtout_sigE_up =   l.diphotonMVA(-1,ipho1,ipho2,ivtx,vtxP,lead_p4,sublead_p4,newSigmaMsUp.first, newSigmaMsUp.second, newSigmaMsUp.first, trainPhil.c_str(), bdtType.c_str(), lead_id_mva, sublead_id_mva);
        float bdtout_sigE_down = l.diphotonMVA(-1,ipho1,ipho2,ivtx,vtxP,lead_p4,sublead_p4,newSigmaMsDown.first, newSigmaMsDown.second, newSigmaMsDown.first, trainPhil.c_str(), bdtType.c_str(), lead_id_mva, sublead_id_mva);
        
        l.FillTree("bdtout",bdtout,"datacard_trees");
        l.FillTree("bdtout_id_up",bdtout_id_up,"datacard_trees");
        l.FillTree("bdtout_id_down",bdtout_id_down,"datacard_trees");
        l.FillTree("bdtout_sigE_up",bdtout_sigE_up,"datacard_trees");
        l.FillTree("bdtout_sigE_down",bdtout_sigE_down,"datacard_trees");
    }
    int proc_id=-1;
    if (proc==Form("ggh_mass_m%3.0f",l.normalizer()->GetMass(cur_type))) proc_id=0;
    if (proc==Form("vbf_mass_m%3.0f",l.normalizer()->GetMass(cur_type))) proc_id=1;
    if (proc==Form("wh_mass_m%3.0f",l.normalizer()->GetMass(cur_type))) proc_id=2;
    if (proc==Form("zh_mass_m%3.0f",l.normalizer()->GetMass(cur_type))) proc_id=3;
    if (proc==Form("tth_mass_m%3.0f",l.normalizer()->GetMass(cur_type))) proc_id=4;
    if (proc==Form("wzh_mass_m%3.0f",l.normalizer()->GetMass(cur_type))) proc_id=5;
    
    l.FillTree("category",category,"datacard_trees");
    l.FillTree("inc_cat",inc_cat,"datacard_trees");
    l.FillTree("process_id",proc_id,"datacard_trees");
    l.FillTree("weight",evweight,"datacard_trees");
    l.FillTree("scale1",scale1,"datacard_trees");
    l.FillTree("scale1_err",scale1_err,"datacard_trees");
    l.FillTree("scale2",scale2,"datacard_trees");
    l.FillTree("scale2_err",scale2_err,"datacard_trees");
    l.FillTree("scale_err",scale_err,"datacard_trees");
    l.FillTree("w_scale_err_2",evweight*scale_err*scale_err,"datacard_trees");
    l.FillTree("smear_err",smear_err,"datacard_trees");
    l.FillTree("w_smear_err_2",evweight*smear_err*smear_err,"datacard_trees");
    l.FillTree("smear1",smear1,"datacard_trees");
    l.FillTree("smear1_err",smear1_err,"datacard_trees");
    l.FillTree("smear2",smear2,"datacard_trees");
    l.FillTree("smear2_err",smear2_err,"datacard_trees");
    l.FillTree("lead_eta",lead_p4.Eta(),"datacard_trees");
    l.FillTree("sublead_eta",sublead_p4.Eta(),"datacard_trees");
    l.FillTree("lead_r9",l.pho_r9[ipho1],"datacard_trees");
    l.FillTree("sublead_r9",l.pho_r9[ipho2],"datacard_trees");
    double lead_calo_eta = ((TVector3*)l.sc_xyz->At(l.pho_scind[ipho1]))->Eta();
    double sublead_calo_eta = ((TVector3*)l.sc_xyz->At(l.pho_scind[ipho2]))->Eta();
    l.FillTree("lead_calo_eta",lead_calo_eta,"datacard_trees");
    l.FillTree("sublead_calo_eta",sublead_calo_eta,"datacard_trees");
    l.FillTree("lead_isEB",TMath::Abs(lead_calo_eta)<1.444,"datacard_trees");
    l.FillTree("sublead_isEB",TMath::Abs(sublead_calo_eta)<1.444,"datacard_trees");
    l.FillTree("vbfmva",myVBF_MVA,"datacard_trees");
    
}

// for Cut-Based
void PhotonAnalysis::saveSpinTree(LoopAll &l, int category, float evweight, TLorentzVector Higgs, TLorentzVector lead_p4, TLorentzVector sublead_p4, int ipho1, int ipho2, int diphoton_id, float vtxProb, bool isCorrectVertex){

   int cur_type = l.itype[l.current];
   TLorentzVector genpho1(0.,0.,0.,0.);
   TLorentzVector genpho2(0.,0.,0.,0.);
   TLorentzVector higgs(0.,0.,0.,0.);

   if (cur_type<0){
        genpho1 = *((TLorentzVector*)l.gh_pho1_p4->At(0));
        genpho2 = *((TLorentzVector*)l.gh_pho2_p4->At(0));
        higgs = *((TLorentzVector*)l.gh_higgs_p4->At(0));
   }
   
   l.FillTree("category",category,"spin_trees");
   l.FillTree("evweight",evweight,"spin_trees");
   
   l.FillTree("higgs_px",Higgs.Px(),"spin_trees");
   l.FillTree("higgs_py",Higgs.Py(),"spin_trees");
   l.FillTree("higgs_pz",Higgs.Pz(),"spin_trees");
   l.FillTree("higgs_E",Higgs.E(),"spin_trees");

   l.FillTree("gh_higgs_px",higgs.Px(),"spin_trees");
   l.FillTree("gh_higgs_py",higgs.Py(),"spin_trees");
   l.FillTree("gh_higgs_pz",higgs.Pz(),"spin_trees");
   l.FillTree("gh_higgs_E",higgs.E(),"spin_trees");

   l.FillTree("lead_px",lead_p4.Px(),"spin_trees");
   l.FillTree("lead_py",lead_p4.Py(),"spin_trees");
   l.FillTree("lead_pz",lead_p4.Pz(),"spin_trees");
   l.FillTree("lead_E",lead_p4.E(),"spin_trees");

   l.FillTree("sublead_px",sublead_p4.Px(),"spin_trees");
   l.FillTree("sublead_py",sublead_p4.Py(),"spin_trees");
   l.FillTree("sublead_pz",sublead_p4.Pz(),"spin_trees");
   l.FillTree("sublead_E",sublead_p4.E(),"spin_trees");

   l.FillTree("gp_lead_px",genpho1.Px(),"spin_trees");
   l.FillTree("gp_lead_py",genpho1.Py(),"spin_trees");
   l.FillTree("gp_lead_pz",genpho1.Pz(),"spin_trees");
   l.FillTree("gp_lead_E",genpho1.E(),"spin_trees");

   l.FillTree("gp_sublead_px",genpho2.Px(),"spin_trees");
   l.FillTree("gp_sublead_py",genpho2.Py(),"spin_trees");
   l.FillTree("gp_sublead_pz",genpho2.Pz(),"spin_trees");
   l.FillTree("gp_sublead_E",genpho2.E(),"spin_trees");

   l.FillTree("costheta_cs",getCosThetaCS(lead_p4,sublead_p4,l.sqrtS),"spin_trees");
   l.FillTree("costheta_hx",getCosThetaHX(lead_p4,sublead_p4,l.sqrtS),"spin_trees");

   l.FillTree("gh_costheta_cs",getCosThetaCS(genpho1,genpho2,l.sqrtS),"spin_trees");
   l.FillTree("gh_costheta_hx",getCosThetaHX(genpho1,genpho2,l.sqrtS),"spin_trees");

   l.FillTree("lead_calo_eta",photonInfoCollection[ipho1].caloPosition().Eta(),"spin_trees");
   l.FillTree("lead_calo_phi",photonInfoCollection[ipho1].caloPosition().Phi(),"spin_trees");
   l.FillTree("lead_r9",l.pho_r9[ipho1],"spin_trees");
   
   l.FillTree("sublead_calo_eta",photonInfoCollection[ipho2].caloPosition().Eta(),"spin_trees");
   l.FillTree("sublead_calo_phi",photonInfoCollection[ipho2].caloPosition().Phi(),"spin_trees");
   l.FillTree("sublead_r9",l.pho_r9[ipho2],"spin_trees");
    
   l.FillTree("rv",isCorrectVertex,"spin_trees");
   l.FillTree("higgs_mass",Higgs.M(),"spin_trees");

    l.FillTree("myVBF_leadEta",myVBFLeadJEta,"spin_trees");
    l.FillTree("myVBF_subleadEta",myVBFSubJEta,"spin_trees");
    l.FillTree("myVBFLeadJPt",myVBFLeadJPt,"spin_trees");
    l.FillTree("myVBFSubJPt",myVBFSubJPt,"spin_trees");
    l.FillTree("myVBF_Mjj",myVBF_Mjj,"spin_trees");

    /*
    PhotonReducedInfo pho1 (
        *((TVector3*)     l.sc_xyz->At(l.pho_scind[ipho1])),
        ((TLorentzVector*)l.pho_p4->At(ipho1))->Energy(),
        energyCorrected[ipho1],
        l.pho_isEB[ipho1], l.pho_r9[ipho1],
        true, // WARNING  setting pass photon ID flag for all photons. This is safe as long as only selected photons are used
        energyCorrectedError[ipho1] // will be altered below, needs to be initialized
    );
    PhotonReducedInfo pho2 (
        *((TVector3*)     l.sc_xyz->At(l.pho_scind[ipho2])),
        ((TLorentzVector*)l.pho_p4->At(ipho2))->Energy(),
        energyCorrected[ipho2],
        l.pho_isEB[ipho2], l.pho_r9[ipho2],
        true, // WARNING  setting pass photon ID flag for all photons. This is safe as long as only selected photons are used
        energyCorrectedError[ipho2] // will be altered below, needs to be initialized
    );
    MassResolution *tempMassRes = new MassResolution;
    tempMassRes->Setup(l,&pho1,&pho2,l.dipho_vtxind[diphoton_id],massResoPars, nR9Categories,nEtaCategories,beamspotSigma,true);

    double sigmaMrv = tempMassRes->massResolutionEonly();
    double sigmaMwv = tempMassRes->massResolutionWrongVtx();
    double lead_sigmaE = tempMassRes->leadPhotonResolution();
    double lead_sigmaE_nosmear = tempMassRes->leadPhotonResolutionNoSmear();
    double sublead_sigmaE = tempMassRes->subleadPhotonResolution();
    double sublead_sigmaE_nosmear = tempMassRes->subleadPhotonResolutionNoSmear();
    
    delete tempMassRes;

   l.FillTree("sigmaMrv",sigmaMrv,"spin_trees");
   l.FillTree("sigmaMwv",sigmaMwv,"spin_trees");
   l.FillTree("lead_sigmaE",lead_sigmaE,"spin_trees");
   l.FillTree("lead_sigmaE_nosmear",lead_sigmaE_nosmear,"spin_trees");
   l.FillTree("sublead_sigmaE",sublead_sigmaE,"spin_trees");
   l.FillTree("sublead_sigmaE_nosmear",sublead_sigmaE_nosmear,"spin_trees");

   l.FillTree("sigmaMoMrv",float(sigmaMrv/Higgs.M()),"spin_trees");
   l.FillTree("sigmaMoMwv",float(sigmaMwv/Higgs.M()),"spin_trees");
   l.FillTree("vtx_prob",vtxProb,"spin_trees");
   l.FillTree("leadPtoM",float(lead_p4.Pt()/Higgs.M()),"spin_trees");
   l.FillTree("subleadPtoM",float(sublead_p4.Pt()/Higgs.M()),"spin_trees");
   l.FillTree("leadEta",float(lead_p4.Eta()),"spin_trees");
   l.FillTree("subleadEta",float(sublead_p4.Eta()),"spin_trees");
   l.FillTree("cosDphi",float(TMath::Cos(lead_p4.Phi()-sublead_p4.Phi())),"spin_trees");
   */
}

// for Mass-factorized
void PhotonAnalysis::saveSpinTree(LoopAll& l, int category, float evweight, TLorentzVector Higgs, TLorentzVector lead_p4, TLorentzVector sublead_p4, int ipho1, int ipho2, float diphobdt, double sigmaMrv, double sigmaMwv, double lead_sigmaE, double lead_sigmaE_nosmear, double sublead_sigmaE, double sublead_sigmaE_nosmear, float vtxProb, float lead_id_mva, float sublead_id_mva){
   
   l.FillTree("category",category,"spin_trees");
   l.FillTree("evweight",evweight,"spin_trees");
   
   l.FillTree("higgs_px",Higgs.Px(),"spin_trees");
   l.FillTree("higgs_py",Higgs.Py(),"spin_trees");
   l.FillTree("higgs_pz",Higgs.Pz(),"spin_trees");
   l.FillTree("higgs_E",Higgs.E(),"spin_trees");

   l.FillTree("lead_px",lead_p4.Px(),"spin_trees");
   l.FillTree("lead_py",lead_p4.Py(),"spin_trees");
   l.FillTree("lead_pz",lead_p4.Pz(),"spin_trees");
   l.FillTree("lead_E",lead_p4.E(),"spin_trees");

   l.FillTree("sublead_px",sublead_p4.Px(),"spin_trees");
   l.FillTree("sublead_py",sublead_p4.Py(),"spin_trees");
   l.FillTree("sublead_pz",sublead_p4.Pz(),"spin_trees");
   l.FillTree("sublead_E",sublead_p4.E(),"spin_trees");

   l.FillTree("costheta_cs",getCosThetaCS(lead_p4,sublead_p4,l.sqrtS),"spin_trees");
   l.FillTree("costheta_hx",getCosThetaHX(lead_p4,sublead_p4,l.sqrtS),"spin_trees");

   l.FillTree("lead_calo_eta",photonInfoCollection[ipho1].caloPosition().Eta(),"spin_trees");
   l.FillTree("lead_calo_phi",photonInfoCollection[ipho1].caloPosition().Phi(),"spin_trees");
   l.FillTree("lead_r9",l.pho_r9[ipho1],"spin_trees");
   
   l.FillTree("sublead_calo_eta",photonInfoCollection[ipho2].caloPosition().Eta(),"spin_trees");
   l.FillTree("sublead_calo_phi",photonInfoCollection[ipho2].caloPosition().Phi(),"spin_trees");
   l.FillTree("sublead_r9",l.pho_r9[ipho2],"spin_trees");

   l.FillTree("diphoton_bdt",diphobdt,"spin_trees");
   l.FillTree("higgs_mass",Higgs.M(),"spin_trees");
 
   l.FillTree("sigmaMrv",sigmaMrv,"spin_trees");
   l.FillTree("sigmaMwv",sigmaMwv,"spin_trees");
   l.FillTree("lead_sigmaE",lead_sigmaE,"spin_trees");
   l.FillTree("lead_sigmaE_nosmear",lead_sigmaE_nosmear,"spin_trees");
   l.FillTree("sublead_sigmaE",sublead_sigmaE,"spin_trees");
   l.FillTree("sublead_sigmaE_nosmear",sublead_sigmaE_nosmear,"spin_trees");

   l.FillTree("sigmaMoMrv",float(sigmaMrv/Higgs.M()),"spin_trees");
   l.FillTree("sigmaMoMwv",float(sigmaMwv/Higgs.M()),"spin_trees");
   l.FillTree("vtx_prob",vtxProb,"spin_trees");
   l.FillTree("leadPtoM",float(lead_p4.Pt()/Higgs.M()),"spin_trees");
   l.FillTree("subleadPtoM",float(sublead_p4.Pt()/Higgs.M()),"spin_trees");
   l.FillTree("leadEta",float(lead_p4.Eta()),"spin_trees");
   l.FillTree("subleadEta",float(sublead_p4.Eta()),"spin_trees");
   l.FillTree("cosDphi",float(TMath::Cos(lead_p4.Phi()-sublead_p4.Phi())),"spin_trees");
   l.FillTree("lead_id_mva",lead_id_mva,"spin_trees");
   l.FillTree("sublead_id_mva",sublead_id_mva,"spin_trees");

}



void PhotonAnalysis::saveVBFTree(LoopAll &l, int category, float evweight, float diphobdt_output){

  l.FillTree("weight",evweight,"vbf_trees");
  l.FillTree("category",category,"vbf_trees");
  l.FillTree("diphomva",diphobdt_output,"vbf_trees");
  l.FillTree("vbfmva",myVBF_MVA,"vbf_trees");
  l.FillTree("leadJPt",myVBFLeadJPt,"vbf_trees");
  l.FillTree("subleadJPt",myVBFSubJPt,"vbf_trees");
  l.FillTree("leadJEta",myVBFLeadJEta,"vbf_trees");
  l.FillTree("subleadJEta",myVBFSubJEta,"vbf_trees");
  l.FillTree("MJJ",myVBF_Mjj,"vbf_trees");
  l.FillTree("dEtaJJ",myVBFdEta,"vbf_trees");
  l.FillTree("zepp",myVBFZep,"vbf_trees");
  l.FillTree("dPhiJJGammaGamma",myVBFdPhi,"vbf_trees");
  l.FillTree("mass",myVBF_Mgg,"vbf_trees");
  
  // gen variables:
  //gen photons,  higgs
  TLorentzVector *genpho1 = (TLorentzVector*)l.gh_pho1_p4->At(0);
  TLorentzVector *genpho2 = (TLorentzVector*)l.gh_pho2_p4->At(0);
  TLorentzVector higgs = *((TLorentzVector*)l.gh_higgs_p4->At(0));
  
  // gen jets
  std::vector<int> sorted_jets;
  for(int ijet=0; ijet<l.genjet_algo1_n; ++ijet) {
      TLorentzVector * p4 = (TLorentzVector*)l.genjet_algo1_p4->At(ijet);
      if( p4->DeltaR( *genpho1 ) > 0.5 && p4->DeltaR( *genpho2 ) > 0.5  ) {
	  sorted_jets.push_back(ijet);
      }
  }
  std::sort(sorted_jets.begin(),sorted_jets.end(),
	    ClonesSorter<TLorentzVector,double,std::greater<double> >(l.genjet_algo1_p4,&TLorentzVector::Pt));
  
  TLorentzVector* j1 ;
  TLorentzVector* j2 ;
  float myVBFdPhi_gen = -99;
  if ( sorted_jets.size() > 1){
      j1 = (TLorentzVector*)l.genjet_algo1_p4->At(sorted_jets[0]);
      j2 = (TLorentzVector*)l.genjet_algo1_p4->At(sorted_jets[1]);
      TLorentzVector dijet = (*j1) + (*j2);
      myVBFdPhi_gen =  fabs(higgs.DeltaPhi(dijet));
  }
  
  
  l.FillTree("dPhiJJGammaGammaGen",myVBFdPhi_gen,"vbf_trees");
}

void PhotonAnalysis::VBFAngles(TLorentzVector& gamma1, TLorentzVector& gamma2, TLorentzVector& J1, TLorentzVector& J2)
{
  myVBFSpin_DeltaPhiJJ = J1.DeltaPhi(J2);
  myVBFSpin_absDeltaPhiJJ = TMath::Abs(myVBFSpin_DeltaPhiJJ);

  std::pair<TLorentzVector, TLorentzVector> IntermediateBoson = GetVBF_IntermediateBoson(gamma1, gamma2, J1, J2);
  TLorentzVector S = IntermediateBoson.first;  //Small deflection
  TLorentzVector L = IntermediateBoson.second; //Large deflection
  TLorentzVector diphoton = gamma1 + gamma2;

  TVector3 boost = -diphoton.BoostVector(); //Get boost vector

  TLorentzVector J1Boosted = J1, J2Boosted = J2, leadingBoosted = gamma1, subleadingBoosted = gamma2, SBoosted = S, LBoosted = L;

  J1Boosted.Boost(boost);
  J2Boosted.Boost(boost);
  leadingBoosted.Boost(boost);
  subleadingBoosted.Boost(boost);
  SBoosted.Boost(boost);
  LBoosted.Boost(boost);

  myVBFSpin_CosThetaJ1 = CosAngle(leadingBoosted, J1Boosted);//leadingBoosted.Vect().Dot(J1Boosted.Vect())/(leadingBoosted.Vect().Mag()*J1Boosted.Vect().Mag());
  myVBFSpin_absCosThetaJ1 = TMath::Abs(myVBFSpin_CosThetaJ1);
  myVBFSpin_CosThetaJ2 = CosAngle(leadingBoosted, J2Boosted);//leadingBoosted.Vect().Dot(J2Boosted.Vect())/(leadingBoosted.Vect().Mag()*J2Boosted.Vect().Mag());
  myVBFSpin_absCosThetaJ2 = TMath::Abs(myVBFSpin_CosThetaJ2);

  // Small deflection
  myVBFSpin_DeltaPhiJJS = GetPerpendicularAngle(SBoosted, J1Boosted, J2Boosted);
  myVBFSpin_absDeltaPhiJJS = TMath::Abs(myVBFSpin_DeltaPhiJJS);
  myVBFSpin_CosThetaS = CosAngle(leadingBoosted, SBoosted);
  myVBFSpin_absCosThetaS = TMath::Abs(myVBFSpin_CosThetaS);
  // Large deflection
  myVBFSpin_DeltaPhiJJL = GetPerpendicularAngle(LBoosted, J1Boosted, J2Boosted);
  myVBFSpin_absDeltaPhiJJL = TMath::Abs(myVBFSpin_DeltaPhiJJL);
  myVBFSpin_CosThetaL = CosAngle(leadingBoosted, LBoosted);
  myVBFSpin_absCosThetaL = TMath::Abs(myVBFSpin_CosThetaL);
}

//Return intermediate boson TLorentzVector for VBF (2 different possible pairings are given
std::pair<TLorentzVector, TLorentzVector> PhotonAnalysis::GetVBF_IntermediateBoson(TLorentzVector& Pho1, TLorentzVector& Pho2, TLorentzVector& Jet1, TLorentzVector& Jet2)
{
  std::pair<TLorentzVector, TLorentzVector> retVal;

  TLorentzVector System = Pho1 + Pho2 + Jet1 + Jet2;
  //std::cout << "M: " << Syst.M() << "; P: " << Syst.P() << "; Pz: " << Syst.Pz() << std::endl;
  Double_t Pz = System.Pz();
  myVBF_Pz = Pz;
  Double_t Mr = System.M2();
  myVBF_S  = Mr;

  //K1 and K2 stand for the incoming parton particles
  Double_t K1_val = (TMath::Sqrt(Pz*Pz+Mr)+Pz)/2.;
  Double_t K2_val = (TMath::Sqrt(Pz*Pz+Mr)-Pz)/2.;
  myVBF_K1 = K1_val;
  myVBF_K2 = K2_val;
  TLorentzVector K1(0, 0,  K1_val, K1_val);
  TLorentzVector K2(0, 0, -K2_val, K2_val);


  std::pair<TLorentzVector*, TLorentzVector*> pairing[2];
  pairing[0].second = &K1;
  pairing[1].second = &K2;
  if(Jet1.Eta() < Jet2.Eta())
  {
    pairing[0].first = &Jet2;
    pairing[1].first = &Jet1;
  }
  else
  {
    pairing[0].first = &Jet1;
    pairing[1].first = &Jet2;
  }

  Int_t index = 1;
  if(TMath::Abs(pairing[0].first->Eta()) < TMath::Abs(pairing[1].first->Eta()))
    index = 0;

  //return value: (S, L)
  //S and L hold the return values (only the most central jet is used):
  // S - pairing using the smallest angle between incoming parton vector and outgoing jet vector
  // L - pairing using the largest angle between incoming parton vector and outgoing jet vector
  retVal.first  = *(pairing[index].second)       - *(pairing[index].first);
  retVal.second = *(pairing[(index+1)%2].second) - *(pairing[index].first);

  return retVal;
}

//Return angle between to TLorentzVectors in the plane perpendicular to a third TLorentzVector
Double_t PhotonAnalysis::GetPerpendicularAngle(TLorentzVector& ref, TLorentzVector& v1, TLorentzVector& v2)
{
  TVector3 ref_norm = (ref.Vect()).Unit();

  TVector3 v1_perp = v1.Vect();
  v1_perp = v1_perp - ref_norm * (ref_norm * v1_perp);

  TVector3 v2_perp = v2.Vect();
  v2_perp = v2_perp - ref_norm * (ref_norm * v2_perp);

  return v1_perp.Angle(v2_perp);
}

double PhotonAnalysis::getCosThetaCS(TLorentzVector g1, TLorentzVector g2,int sqrtS){
    
    TLorentzVector b1,b2,diphoton;
    b1.SetPx(0); b1.SetPy(0);
    b2.SetPx(0); b2.SetPy(0);
    double beamE = 500.*sqrtS;
    b1.SetPz( beamE); b1.SetE(beamE);
    b2.SetPz(-beamE); b2.SetE(beamE);
    
    diphoton=g1+g2;
    TVector3 boostToDiphotonFrame = -diphoton.BoostVector();

    // Boost to higgs frame
    TLorentzVector refDIPHO_g1 = g1; refDIPHO_g1.Boost(boostToDiphotonFrame);
    TLorentzVector refDIPHO_b1 = b1; refDIPHO_b1.Boost(boostToDiphotonFrame);
    TLorentzVector refDIPHO_b2 = b2; refDIPHO_b2.Boost(boostToDiphotonFrame);

    // Getting beam 3-vector from 4-vectors
    TVector3 refDIPHO_vb1_direction = refDIPHO_b1.Vect().Unit();
    TVector3 refDIPHO_vb2_direction = refDIPHO_b2.Vect().Unit();

    // Definition of zz directions
    TVector3 direction_cs = (refDIPHO_vb1_direction - refDIPHO_vb2_direction).Unit(); // CS direction

    return TMath::Cos(direction_cs.Angle(refDIPHO_g1.Vect()));
}


float PhotonAnalysis::getDiphoBDTOutput(LoopAll &l,int diphoton_id, TLorentzVector lead_p4,TLorentzVector sublead_p4, std::string bdtTrainingPhilosophy){

    //    cout<<"[DEBUG]: index"<<l.dipho_leadind[diphoton_id]<<" index sublead "<<l.dipho_subleadind[diphoton_id]<<" nphot "<<l.pho_n<<endl;
    //      cout<<"[DEBUG]:pt"<<lead_p4.Pt()<<endl;                                                                                                                     
    
    //setting resolutions
    massResolutionCalculator->Setup(l,&photonInfoCollection[l.dipho_leadind[diphoton_id]],&photonInfoCollection[l.dipho_subleadind[diphoton_id]],0,//default vertex
				    massResoPars,nR9Categories,nEtaCategories,beamspotSigma,true);

    float sigmaMrv = massResolutionCalculator->massResolutionCorrVtx();
    float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
    float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
    
    //diphoton mva                                                                                                                                                     
    float diphobdt_output = l.diphotonMVA(-1,l.dipho_leadind[diphoton_id],l.dipho_subleadind[diphoton_id],0 ,//vertex 0 probability 1                             
                                          1,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,
                                          bdtTrainingPhilosophy.c_str(), bdtTrainingType.c_str(),
                                          -1.,-1.);

    return diphobdt_output;

}

double PhotonAnalysis::getCosThetaHX(TLorentzVector g1, TLorentzVector g2,int sqrtS){
    
    TLorentzVector b1,b2,diphoton;
    b1.SetPx(0); b1.SetPy(0);
    b2.SetPx(0); b2.SetPy(0);
    double beamE = 500.*sqrtS;
    b1.SetPz( beamE); b1.SetE(beamE);
    b2.SetPz(-beamE); b2.SetE(beamE);

  diphoton=g1+g2;
  TVector3 boostToDiphotonFrame = -diphoton.BoostVector();

  // Boost to higgs frame
  TLorentzVector refDIPHO_g1 = g1; refDIPHO_g1.Boost(boostToDiphotonFrame);

  // Definition of zz directions
  TVector3 direction_hx = diphoton.Vect().Unit(); // HX direction

  return TMath::Cos(direction_hx.Angle(refDIPHO_g1.Vect()));

}

void PhotonAnalysis::VHLepTag2013(LoopAll& l, int & diphotonVHlep_id, bool & VHlep1event, bool & VHlep2event, bool mvaselection, int & mu_ind, int & muVtx, int VHmuevent_cat, int & el_ind, int & elVtx, int VHelevent_cat, float* smeared_pho_energy, float phoidMvaCut, float eventweight, std::vector<float> smeared_pho_weight, bool isSyst, bool vetodipho, bool kinonly){
    bool VHmuevent_prov=false;
    bool VHelevent_prov=false;
    if(mvaselection){
        VHmuevent_prov=MuonTag2012B(l,diphotonVHlep_id,mu_ind,muVtx,VHmuevent_cat,&smeared_pho_energy[0],lep_sync,mvaselection,phoidMvaCut,eventweight,smeared_pho_weight, !isSyst, vetodipho, kinonly);
        int diphotonVH_ele_id=-1;
        VHelevent_prov=ElectronTag2012B(l,diphotonVH_ele_id,el_ind,elVtx,VHelevent_cat,&smeared_pho_energy[0],lep_sync,mvaselection,phoidMvaCut,eventweight,smeared_pho_weight, !isSyst, vetodipho, kinonly);
        if(!VHmuevent_prov && VHelevent_prov) diphotonVHlep_id=diphotonVH_ele_id;
    } else {
        VHmuevent_prov=MuonTag2012B(l,diphotonVHlep_id,mu_ind,muVtx,VHmuevent_cat,&smeared_pho_energy[0],lep_sync,false,-0.2,eventweight,smeared_pho_weight,!isSyst, vetodipho, kinonly);
        int diphotonVH_ele_id=-1;
        VHelevent_prov=ElectronTag2012B(l,diphotonVH_ele_id,el_ind,elVtx,VHelevent_cat,&smeared_pho_energy[0],lep_sync,false,-0.2,eventweight,smeared_pho_weight,!isSyst, vetodipho, kinonly);
        if(!VHmuevent_prov && VHelevent_prov) diphotonVHlep_id=diphotonVH_ele_id;
    }
    int vertex = l.dipho_vtxind[diphotonVHlep_id];
    if(VHmuevent_prov || VHelevent_prov){
        int Njet_lepcat = VHNumberOfJets(l, diphotonVHlep_id, vertex, VHelevent_prov, VHmuevent_prov, el_ind, mu_ind, &smeared_pho_energy[0]);
        if(Njet_lepcat<3) l.VHNewLeptonCategorization(VHlep1event, VHlep2event, diphotonVHlep_id, vertex, VHelevent_prov, VHmuevent_prov, el_ind, mu_ind, &smeared_pho_energy[0], 45.0, moriond2013MetCorrection);
    }
    l.VHTwoMuonsEvents(VHlep1event, VHlep2event, diphotonVHlep_id, muVtx, &smeared_pho_energy[0], leadEtVHlepCut, subleadEtVHlepCut, applyPtoverM, mvaselection, diphobdt_output_Cut_VHLep, phoidMvaCut, vetodipho, kinonly, bdtTrainingType.c_str());
    l.VHTwoElectronsEvents(VHlep1event, VHlep2event, diphotonVHlep_id, elVtx, &smeared_pho_energy[0], leadEtVHlepCut, subleadEtVHlepCut, applyPtoverM, mvaselection, diphobdt_output_Cut_VHLep, phoidMvaCut, vetodipho, kinonly, bdtTrainingType.c_str());
}

int PhotonAnalysis::VHNumberOfJets(LoopAll& l, int diphotonVHlep_id, int vertex, bool VHelevent_prov, bool VHmuevent_prov, int el_ind, int mu_ind, float* smeared_pho_energy){

  TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], vertex, &smeared_pho_energy[0]);
  TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], vertex, &smeared_pho_energy[0]);  
  int Njet_lepcat = 0;

  bool * jetid_flags=0; //PU-JET VETO
  static std::vector<unsigned char> id_flags;
  switchJetIdVertex( l, l.dipho_vtxind[diphotonVHlep_id] );
  id_flags.resize(l.jet_algoPF1_n);
  for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
    id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
  }
  jetid_flags = (bool*)&id_flags[0];

  for(int i=0; i<l.jet_algoPF1_n; i++){
    TLorentzVector * p4_jet = (TLorentzVector *) l.jet_algoPF1_p4->At(i);
    double dR_jet_PhoLead = p4_jet->DeltaR(lead_p4);
    double dR_jet_PhoSubLead = p4_jet->DeltaR(sublead_p4);
    double dR_jet_muon = 10.0;
    double dR_jet_electron = 10.0;
    
    if(VHelevent_prov){ 
      TLorentzVector* el_jet = (TLorentzVector*) l.el_std_p4->At(el_ind); 
      dR_jet_electron = p4_jet->DeltaR(*el_jet);
    }
    if(VHmuevent_prov){ 
      TLorentzVector* mu_jet = (TLorentzVector*) l.mu_glo_p4->At(mu_ind); 
      dR_jet_muon = p4_jet->DeltaR(*mu_jet);
    }
    if(dR_jet_PhoLead<0.5) continue;
    if(dR_jet_PhoSubLead<0.5) continue;
    if(dR_jet_electron<0.5) continue;
    if(dR_jet_muon<0.5) continue;
    if(p4_jet->Eta()>2.4) continue;
    if(p4_jet->Pt()<20) continue;
    if(jetid_flags != 0 && !jetid_flags[i]) continue;  //PILEUP
    Njet_lepcat = Njet_lepcat + 1;
  }

  return Njet_lepcat;
}

void PhotonAnalysis::GetRegressionCorrections(LoopAll &l){

    if (regressionVersion==5){
        GetRegressionCorrectionsV5(l);
    } else if (regressionVersion==8) {
        // v6/v7 used for 7 TeV regression (v6 Barrel, v7 Endcap)
        // handled in V8
        GetRegressionCorrectionsV8(l);      
    }
}
void PhotonAnalysis::GetRegressionCorrectionsV8(LoopAll &l){
    // V7 7TeV Endcap use
    for (int ipho=0;ipho<l.pho_n;ipho++){
        double ecor,ecorerr;
    
        int sc_index = l.pho_scind[ipho];

        TVector3 *sc = ((TVector3*)l.sc_xyz->At(sc_index)); 
        bool isbarrel = (fabs(sc->Eta())<1.48);
     
        if (isbarrel){
          GetSinglePhotonRegressionCorrectionV6(l,ipho,&ecor,&ecorerr);      
        } else {
          GetSinglePhotonRegressionCorrectionV7(l,ipho,&ecor,&ecorerr);      
        }
        // Save new branches 
        l.pho_regr_energy_otf[ipho] = ecor;
        l.pho_regr_energyerr_otf[ipho] = ecorerr;
    }
}
void PhotonAnalysis::GetSinglePhotonRegressionCorrectionV7(LoopAll &l, int ipho, double *ecor, double *ecorerr){
    // V7 7TeV Endcap use

    double cbalpha1,cbn1,cbalpha2,cbn2,pdfpeakval,cbmean,cbsigma;

    double phoE = ((TLorentzVector*)l.pho_p4->At(ipho))->Energy();
    double r9=l.pho_r9[ipho];

    int sc_index      = l.pho_scind[ipho];
    TVector3 *sc = ((TVector3*)l.sc_xyz->At(sc_index)); 

    int sc_seed_index = l.sc_bcseedind[sc_index];

    TVector3 *bcpos =(TVector3*)l.bc_xyz->At(sc_seed_index);
    double bcE = ((TLorentzVector*)l.bc_p4->At(sc_seed_index))->Energy();

    //   //basic supercluster variables
    _vals[0]  = l.sc_raw[sc_index];
    _vals[1]  = sc->Eta();
    _vals[2]  = r9;
    _vals[3] = l.sc_seta[sc_index];
    _vals[4] = l.sc_sphi[sc_index];
    _vals[5] = (double)l.sc_nbc[sc_index];
    _vals[6] = l.pho_hoe_bc[ipho];//p.hadTowOverEm();
    _vals[7] = l.rho_algo1;
    _vals[8] = (double)l.vtx_std_n;//double(vtxcol.size());

    //seed basic cluster variables
    double bemax = l.bc_s1[sc_seed_index];//clustertools.eMax(*b);
    double be2nd = l.pho_e2nd[ipho];//clustertools.e2nd(*b);
    double betop = l.pho_etop[ipho];//clustertools.eTop(*b);
    double bebottom = l.pho_ebottom[ipho];//clustertools.eBottom(*b);
    double beleft = l.pho_eleft[ipho];//clustertools.eLeft(*b);
    double beright = l.pho_eright[ipho];//clustertools.eRight(*b);

    double be2x5max = l.pho_e2x5max[ipho];//clustertools.e2x5Max(*b);
    double be2x5top = l.pho_e2x5top[ipho];//clustertools.e2x5Top(*b);
    double be2x5bottom = l.pho_e2x5bottom[ipho];//clustertools.e2x5Bottom(*b);
    double be2x5left = l.pho_e2x5left[ipho];//clustertools.e2x5Left(*b);
    double be2x5right = l.pho_e2x5right[ipho];//clustertools.e2x5Right(*b);

    double be5x5 = l.bc_s25[sc_seed_index];//clustertools.e5x5(*b);
    double be3x3 = l.bc_s9[sc_seed_index];//clustertools.e5x5(*b);

    _vals[9] = bcpos->Eta()-sc->Eta();
    _vals[10] = l.DeltaPhi(bcpos->Phi(),sc->Phi());
    _vals[11] = bcE/l.sc_raw[sc_index];
    _vals[12] = be3x3/be5x5;
    _vals[13] = l.bc_sieie[sc_seed_index]; //sigietaieta (this is stored in bc collection)
    _vals[14] = TMath::Sqrt(l.pho_sipip[ipho]); //sigiphiiphi
    _vals[15] = l.pho_sieip[ipho];//clustertools.localCovariances(*b)[1];       //sigietaiphi

    _vals[16] = bemax/be5x5;                       //crystal energy ratio gap variables   
    _vals[17] = be2nd/be5x5;
    _vals[18] = betop/be5x5;
    _vals[19] = bebottom/be5x5;
    _vals[20] = beleft/be5x5;
    _vals[21] = beright/be5x5;
    _vals[22] = be2x5max/be5x5;                       //crystal energy ratio gap variables   
    _vals[23] = be2x5top/be5x5;
    _vals[24] = be2x5bottom/be5x5;
    _vals[25] = be2x5left/be5x5;
    _vals[26] = be2x5right/be5x5;

    // Always assume endcap for V7
    //preshower energy ratio (endcap only)
    _vals[27]  = l.sc_pre[sc_index]/l.sc_raw[sc_index];

    double den =  l.sc_pre[sc_index]+l.sc_raw[sc_index];
      
    //set raw response variables from GBRForest
    _sigma->setVal(_forestDee->GetResponse(&_vals[0],0));
    _mean->setVal(_forestDee->GetResponse(&_vals[0],1));
    _n1->setVal(_forestDee->GetResponse(&_vals[0],2));
    _n2->setVal(_forestDee->GetResponse(&_vals[0],3));
    
    //retrieve final pdf parameter values from transformed forest outputs
    cbmean = _meanlim->getVal();
    cbsigma = _sigmalim->getVal();
    cbalpha1 = 2.0;  //alpha hardcoded in this version of the regression
    cbn1 = _n1lim->getVal();
    cbalpha2 = 1.0;  //alpha hardcoded in this version of the regression
    cbn2 = _n2lim->getVal();
    
    _tgt->setVal(cbmean); //evaluate pdf at peak position  
    pdfpeakval = _pdf->getVal(*_tgt);
      
    //set final energy and relative energy resolution
    *ecor = den*cbmean;
    double sigEoverE = cbsigma/cbmean;
    *ecorerr = sigEoverE*(*ecor);


    //// // Set vectors used in reduction;
    //// energyCorrected[ipho] = ecor;
    //// energyCorrectedError[ipho] = ecorerr;



    //// // Overwrite old branches
    //// l.pho_regr_energy[ipho] = ecor;
    //// l.pho_regr_energyerr[ipho] = ecorerr;
    
}

void PhotonAnalysis::GetSinglePhotonRegressionCorrectionV6(LoopAll &l, int ipho, double *ecor, double *ecorerr){
   
    // V6 7TeV Barrel use  

    double cbalpha1,cbn1,cbalpha2,cbn2,pdfpeakval,cbmean,cbsigma;

    double phoE = ((TLorentzVector*)l.pho_p4->At(ipho))->Energy();
    double r9=l.pho_r9[ipho];

    int sc_index      = l.pho_scind[ipho];
    TVector3 *sc = ((TVector3*)l.sc_xyz->At(sc_index)); 

    int sc_seed_index = l.sc_bcseedind[sc_index];

    TVector3 *bcpos =(TVector3*)l.bc_xyz->At(sc_seed_index);
    double bcE = ((TLorentzVector*)l.bc_p4->At(sc_seed_index))->Energy();

    //   //basic supercluster variables
    _vals[0]  = l.sc_raw[sc_index];
    _vals[1]  = sc->Eta();
    _vals[2]  = sc->Phi();
    _vals[3]  = r9;
    _vals[4] = l.sc_seta[sc_index];
    _vals[5] = l.sc_sphi[sc_index];
    _vals[6] = (double)l.sc_nbc[sc_index];
    _vals[7] = l.pho_hoe_bc[ipho];//p.hadTowOverEm();
    _vals[8] = l.rho_algo1;
    _vals[9] = (double)l.vtx_std_n;//double(vtxcol.size());

    //seed basic cluster variables
    double bemax = l.bc_s1[sc_seed_index];//clustertools.eMax(*b);
    double be2nd = l.pho_e2nd[ipho];//clustertools.e2nd(*b);
    double betop = l.pho_etop[ipho];//clustertools.eTop(*b);
    double bebottom = l.pho_ebottom[ipho];//clustertools.eBottom(*b);
    double beleft = l.pho_eleft[ipho];//clustertools.eLeft(*b);
    double beright = l.pho_eright[ipho];//clustertools.eRight(*b);

    double be2x5max = l.pho_e2x5max[ipho];//clustertools.e2x5Max(*b);
    double be2x5top = l.pho_e2x5top[ipho];//clustertools.e2x5Top(*b);
    double be2x5bottom = l.pho_e2x5bottom[ipho];//clustertools.e2x5Bottom(*b);
    double be2x5left = l.pho_e2x5left[ipho];//clustertools.e2x5Left(*b);
    double be2x5right = l.pho_e2x5right[ipho];//clustertools.e2x5Right(*b);

    double be5x5 = l.bc_s25[sc_seed_index];//clustertools.e5x5(*b);
    double be3x3 = l.bc_s9[sc_seed_index];//clustertools.e5x5(*b);

    _vals[10] = bcpos->Eta()-sc->Eta();
    _vals[11] = l.DeltaPhi(bcpos->Phi(),sc->Phi());
    _vals[12] = bcE/l.sc_raw[sc_index];
    _vals[13] = be3x3/be5x5;
    _vals[14] = l.bc_sieie[sc_seed_index]; //sigietaieta (this is stored in bc collection)
    _vals[15] = TMath::Sqrt(l.pho_sipip[ipho]); //sigiphiiphi
    _vals[16] = l.pho_sieip[ipho];//clustertools.localCovariances(*b)[1];       //sigietaiphi

    _vals[17] = bemax/be5x5;                       //crystal energy ratio gap variables   
    _vals[18] = be2nd/be5x5;
    _vals[19] = betop/be5x5;
    _vals[20] = bebottom/be5x5;
    _vals[21] = beleft/be5x5;
    _vals[22] = beright/be5x5;
    _vals[23] = be2x5max/be5x5;                       //crystal energy ratio gap variables   
    _vals[24] = be2x5top/be5x5;
    _vals[25] = be2x5bottom/be5x5;
    _vals[26] = be2x5left/be5x5;
    _vals[27] = be2x5right/be5x5;

    // V6 Is always Barrel so use that implementation !
        //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    _vals[28] = be5x5/bcE;

    int bieta = l.pho_bieta[ipho];
    int biphi = l.pho_biphi[ipho];     

    _vals[29] = bieta; //crystal ieta
    _vals[30] = biphi; //crystal iphi
    _vals[31] = (bieta-1*std::abs(bieta)/bieta)%5;; //submodule boundary eta symmetry
    _vals[32] = (biphi-1)%2; //submodule boundary phi symmetry
    _vals[33] = (TMath::Abs(bieta)<=25)*((bieta-1*TMath::Abs(bieta)/bieta)%25) + (TMath::Abs(bieta)>25)*((bieta-26*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry
    _vals[34] = biphi%20; //module boundary phi symmetry
    _vals[35] = l.pho_betacry[ipho];//betacry; //local coordinates with respect to closest crystal center at nominal shower depth
    _vals[36] = l.pho_phicry[ipho];//bphicry;


    double den = l.sc_raw[sc_index];
    //set raw response variables from GBRForest
    _sigma->setVal(_forestDeb->GetResponse(&_vals[0],0));
    _mean->setVal(_forestDeb->GetResponse(&_vals[0],1));
    _n1->setVal(_forestDeb->GetResponse(&_vals[0],2));
    _n2->setVal(_forestDeb->GetResponse(&_vals[0],3));
  
    //retrieve final pdf parameter values from transformed forest outputs
    cbmean = _meanlim->getVal();
    cbsigma = _sigmalim->getVal();
    cbalpha1 = 2.0;  //alpha hardcoded in this version of the regression
    cbn1 = _n1lim->getVal();
    cbalpha2 = 1.0;  //alpha hardcoded in this version of the regression
    cbn2 = _n2lim->getVal();
  
    _tgt->setVal(cbmean); //evaluate pdf at peak position  
    pdfpeakval = _pdf->getVal(*_tgt);
    
    //set final energy and relative energy resolution
    *ecor = den*cbmean;
    double sigEoverE = cbsigma/cbmean;
    *ecorerr  = sigEoverE*(*ecor); // note difference from V5 convention
    
}

void PhotonAnalysis::GetRegressionCorrectionsV5(LoopAll &l){

    // v5 used for 8 TeV Energy Regression 
    // On the fly energy regression values
    for (int ipho=0;ipho<l.pho_n;ipho++){

        double ecor,ecorerr,cbalpha1,cbn1,cbalpha2,cbn2,pdfpeakval;

        double phoE = ((TLorentzVector*)l.pho_p4->At(ipho))->Energy();
        double r9=l.pho_r9[ipho];

        int sc_index      = l.pho_scind[ipho];
        int sc_seed_index = l.sc_bcseedind[sc_index];

        TVector3 *sc = ((TVector3*)l.sc_xyz->At(sc_index)); 
        TVector3 *bcpos =(TVector3*)l.bc_xyz->At(sc_seed_index);
        double bcE = ((TLorentzVector*)l.bc_p4->At(sc_seed_index))->Energy();


        // New semi-parametric regression 
        bool isbarrel = (fabs(sc->Eta())<1.48); 

        //   //basic supercluster variables
        _vals[0]  = l.sc_raw[sc_index];
        _vals[1]  = sc->Eta();
        _vals[2]  = r9;
        _vals[3] = l.sc_seta[sc_index];
        _vals[4] = l.sc_sphi[sc_index];
        _vals[5] = (double)l.sc_nbc[sc_index];
        _vals[6] = l.pho_hoe_bc[ipho];//p.hadTowOverEm();
        _vals[7] = l.rho_algo1;
        _vals[8] = (double)l.vtx_std_n;//double(vtxcol.size());

        //seed basic cluster variables
        double bemax = l.bc_s1[sc_seed_index];//clustertools.eMax(*b);
        double be2nd = l.pho_e2nd[ipho];//clustertools.e2nd(*b);
        double betop = l.pho_etop[ipho];//clustertools.eTop(*b);
        double bebottom = l.pho_ebottom[ipho];//clustertools.eBottom(*b);
        double beleft = l.pho_eleft[ipho];//clustertools.eLeft(*b);
        double beright = l.pho_eright[ipho];//clustertools.eRight(*b);

        double be2x5max = l.pho_e2x5max[ipho];//clustertools.e2x5Max(*b);
        double be2x5top = l.pho_e2x5top[ipho];//clustertools.e2x5Top(*b);
        double be2x5bottom = l.pho_e2x5bottom[ipho];//clustertools.e2x5Bottom(*b);
        double be2x5left = l.pho_e2x5left[ipho];//clustertools.e2x5Left(*b);
        double be2x5right = l.pho_e2x5right[ipho];//clustertools.e2x5Right(*b);

        double be5x5 = l.bc_s25[sc_seed_index];//clustertools.e5x5(*b);
        double be3x3 = l.bc_s9[sc_seed_index];//clustertools.e5x5(*b);

        _vals[9] = bcpos->Eta()-sc->Eta();
        _vals[10] = bcpos->DeltaPhi(*sc);
        _vals[11] = bcE/l.sc_raw[sc_index];
        _vals[12] = be3x3/be5x5;
        _vals[13] = l.bc_sieie[sc_seed_index]; //sigietaieta (this is stored in bc collection)
        _vals[14] = TMath::Sqrt(l.pho_sipip[ipho]); //sigiphiiphi
        _vals[15] = l.pho_sieip[ipho];//clustertools.localCovariances(*b)[1];       //sigietaiphi

        _vals[16] = bemax/be5x5;                       //crystal energy ratio gap variables   
        _vals[17] = be2nd/be5x5;
        _vals[18] = betop/be5x5;
        _vals[19] = bebottom/be5x5;
        _vals[20] = beleft/be5x5;
        _vals[21] = beright/be5x5;
        _vals[22] = be2x5max/be5x5;                       //crystal energy ratio gap variables   
        _vals[23] = be2x5top/be5x5;
        _vals[24] = be2x5bottom/be5x5;
        _vals[25] = be2x5left/be5x5;
        _vals[26] = be2x5right/be5x5;

        if (isbarrel) {
            //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
            _vals[27] = be5x5/bcE;

            int bieta = l.pho_bieta[ipho];
            int biphi = l.pho_biphi[ipho];     

            _vals[28] = bieta; //crystal ieta
            _vals[29] = biphi%18; //crystal iphi supermodule symmetry
            _vals[30] = bieta%5; //submodule boundary eta symmetry
            _vals[31] = biphi%2; //submodule boundary phi symmetry
            _vals[32] = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry
            _vals[33] = biphi%20; //module boundary phi symmetry
            _vals[34] = l.pho_betacry[ipho];//betacry; //local coordinates with respect to closest crystal center at nominal shower depth
            _vals[35] = l.pho_phicry[ipho];//bphicry;
            
        }
        else {
            //preshower energy ratio (endcap only)
            _vals[27]  = l.sc_pre[sc_index]/l.sc_raw[sc_index];
            
        }

        double den;
        HybridGBRForest *forest;  
        if (isbarrel) {
            den = l.sc_raw[sc_index];
            forest = _foresteb;
        }
        else {
            den = l.sc_raw[sc_index]+l.sc_pre[sc_index];
            forest = _forestee;
        }

        _tgt->setVal(1.0); //evaluate pdf at peak position

        //set raw response variables from GBRForest
        _sigma->setVal(forest->GetResponse(&_vals[0],0));
        _mean->setVal(forest->GetResponse(&_vals[0],1));
        _n1->setVal(forest->GetResponse(&_vals[0],2));
        _n2->setVal(forest->GetResponse(&_vals[0],3));

        //retrieve final pdf parameter values from transformed forest outputs
        // cbsigma (sigmalim->getVal()) is the sigmaE/E so in the branch we save cbsigma*ecor (ie the absolute error in GeV) 
        ecor = den/_meanlim->getVal();
        ecorerr = _sigmalim->getVal()*ecor;

        cbalpha1 = 2.0;  //alpha hardcoded in this version of the regression
        cbn1 = _n1lim->getVal();
        cbalpha2 = 1.0;  //alpha hardcoded in this version of the regression
        cbn2 = _n2lim->getVal();

        // note ecor is now the corrected energy (save this directly)
        pdfpeakval = _pdf->getVal(*_tgt);


        //// // Set vectors used in reduction;
        //// energyCorrected[ipho] = ecor;
        //// energyCorrectedError[ipho] = ecorerr;

        // Save new branches 
        l.pho_regr_energy_otf[ipho] = ecor;
        l.pho_regr_energyerr_otf[ipho] = ecorerr;

        if (PADEBUG) {
            std::cout << "PhotonAnalysis::GetRegressionCorrections ----/" <<std::endl;
            std::cout << " Is barrel? .... " << isbarrel <<std::endl;
            std::cout << " Inputs ...." << std::endl;
            for (int vi=0;vi<36;vi++){ //36 params for Photon correction
                std::cout << "Val " << vi << " " << _vals[vi] << std::endl;
            }
            std::cout << " photon Energy in ntuple / new value     " << ipho << " = " << l.pho_regr_energy[ipho] << " / " << ecor  << std::endl;
            std::cout << " photon Resolution in ntuple / new value " << ipho << " = " << l.pho_regr_energyerr[ipho] << " / " << ecorerr << std::endl; 
            std::cout << "---------------------------------------------/" <<std::endl;
        }

        //// // Overwrite old branches
        //// l.pho_regr_energy[ipho] = ecor;
        //// l.pho_regr_energyerr[ipho] = ecorerr;
    }
}


std::pair<int, int> PhotonAnalysis::SelectBtaggedAndHighestPtJets(LoopAll& l,int diphoton_id, const TLorentzVector& leadpho,const TLorentzVector& subleadpho, Bool_t * jetid_flags)
{
    std::pair<int, int> myJets(-1,-1);
    std::pair<int, int> fail(-1,-1);

    std::pair<float, float> myJetspt(-1.,-1.);

    float dr2pho = 0.5;
    float dr2jet = 0.5;

    TLorentzVector* j1p4;
    TLorentzVector* j2p4;
    float j1pt=-1;
    float j2pt=-1;

    float ptJets_thresh=25.;

    static std::vector<unsigned char> id_flags;
    if( jetid_flags == 0 ) {
      ((PhotonAnalysis*) NULL)->switchJetIdVertex( l, l.dipho_vtxind[diphoton_id] );
      id_flags.resize(l.jet_algoPF1_n);
      for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
	id_flags[ijet] = PileupJetIdentifier::passJetId(l.jet_algoPF1_cutbased_wp_level[ijet], PileupJetIdentifier::kLoose);
      }
      jetid_flags = (bool*)&id_flags[0];
    }



    // select btagged or highest pt jets
    std::vector<int> index_selected_btagloose;
    for(int j1_i=0; j1_i<l.jet_algoPF1_n; j1_i++){
        j1p4 = (TLorentzVector*) l.jet_algoPF1_p4->At(j1_i);
        if(jetid_flags != 0 && !jetid_flags[j1_i]) continue; 
        if(fabs(j1p4->Eta()) > 2.4) continue;
        if(j1p4->DeltaR(leadpho) < dr2pho) continue;
        if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
        j1pt=j1p4->Pt();
	if(j1pt<ptJets_thresh) continue;

	if(l.jet_algoPF1_csvBtag[j1_i]>0.244) {
	  index_selected_btagloose.push_back(j1_i);
	  }
	
	  if(j1pt>myJetspt.first) {
            myJets.second=myJets.first;
            myJetspt.second=myJetspt.first;
            myJetspt.first=j1pt;
            myJets.first=j1_i;
	  }        else if(j1pt>myJetspt.second) {
            myJetspt.second=j1pt;
            myJets.second=j1_i;
        }

    }
    
    if( index_selected_btagloose.size()==1 ) {
      if( index_selected_btagloose[0]!=myJets.first ) {
	myJets.second = index_selected_btagloose[0];

      } 
 
   }else if (index_selected_btagloose.size()>1){
      myJets.first = index_selected_btagloose[0];
      myJets.second = index_selected_btagloose[1];
    }


    return myJets;

}



// Local Variables:
// mode: c++
// c-basic-offset: 4
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
