#include "../interface/PhotonPlusJetVertexAnalysis.h"
#include "PhotonReducedInfo.h"
#define PADEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
PhotonPlusJetVertexAnalysis::PhotonPlusJetVertexAnalysis()  
{
    name_ = "PhotonPlusJetVertexAnalysis";
}

// ----------------------------------------------------------------------------------------------------
PhotonPlusJetVertexAnalysis::~PhotonPlusJetVertexAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void PhotonPlusJetVertexAnalysis::Term(LoopAll& l) 
{

}

// ----------------------------------------------------------------------------------------------------
void PhotonPlusJetVertexAnalysis::Init(LoopAll& l) 
{  
  doSystematics = false;
  StatAnalysis::Init(l);
}


// ----------------------------------------------------------------------------------------------------
bool PhotonPlusJetVertexAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4,
                         float & mass, float & evweight, int & category, int & diphoton_id, bool & isCorrectVertex, float &kinematic_bdtout,
                         bool isSyst, 
                         float syst_shift, bool skipSelection, 
                         BaseGenLevelSmearer *genSys, BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys) 
{
  // check trigger
  bool passTrigger=false;
  for (unsigned int j=0; j<(unsigned int) l.hlt_bit->size(); j++) {
    TString TriggerName(l.hlt_path_names_HLT->at(l.hlt_bit->at(j)));
    if (TriggerName.Contains("HLT_Photon50_CaloIdVL_IsoL_v")) passTrigger=true;
  }
  if (!passTrigger) return false;


  if (l.pho_n < 1) return false;

  int cur_type = l.itype[l.current];
  //float sampleweight = l.sampleContainer[l.current_sample_index].weight;
 
  // first apply corrections and smearing on the single photons 
  if (PADEBUG) cout << "[PADEBUG] : applying energy smearings" << endl;
  smeared_pho_energy.clear(); smeared_pho_energy.resize(l.pho_n,0.); 
  smeared_pho_r9.clear();     smeared_pho_r9.resize(l.pho_n,0.); 
  smeared_pho_weight.clear(); smeared_pho_weight.resize(l.pho_n,1.);
  applySinglePhotonSmearings(smeared_pho_energy, smeared_pho_r9, smeared_pho_weight, cur_type, l, energyCorrected, energyCorrectedError,
			     phoSys, syst_shift);
    
  int vtxind = 0;
  
  // photon selection : 
  // -----  pT > 30 GeV
  // -----  pass CiC supertight
  // if more than one photon, take the one with highest pt (maybe we could use all the photons....)
  if (PADEBUG) cout << "[PADEBUG] : starting photon selection" << endl;  
  float maxphopt  = -1;
  float maxphoind = 0 ;
  for(size_t ipho=0; ipho<l.pho_n; ipho++ ) {
    TLorentzVector phoP4 = l.get_pho_p4( ipho, vtxind, &smeared_pho_energy[0] ); 
    if ( phoP4.Pt() < 30 ) continue;
    if ( (l.pho_cic4pfcutlevel_lead)->at(ipho).at(vtxind)<4) continue;
    if (PADEBUG) cout << "   photon  pt :" << phoP4.Pt() << endl;
    if (PADEBUG) cout << "   photon  id :" << (l.pho_cic4pfcutlevel_lead)->at(ipho).at(0) << endl;
    if ( phoP4.Pt() > maxphopt ) {
      maxphopt = phoP4.Pt();
      maxphoind = ipho;
    }
  }    
 
  if ( maxphoind <0 ) return false;

  if (PADEBUG) cout << "[PADEBUG] : found one good photon" << endl;
  
  TLorentzVector photon = l.get_pho_p4( maxphoind, vtxind, &smeared_pho_energy[0] ); 
  

  // jet selection
  if (PADEBUG) cout << "[PADEBUG] : checking jets..." << endl;
  int maxjetind=-1;
  TLorentzVector MaxJetP4(0,0,0,0);
  TLorentzVector MaxJetTrackP4(0,0,0,0);
  for (unsigned int ijet=0; ijet < (unsigned int)l.jet_algoPF1_n; ijet++) {
    TLorentzVector JetP4 = *((TLorentzVector*)l.jet_algoPF1_p4->At(ijet));
    TLorentzVector TrackSumP4(0,0,0,0);
    double dR = JetP4.DeltaR(photon);
    if (JetP4.Eta()>2.4 || dR<0.4 || l.jet_algoPF1_ntk[ijet]<3) continue;
    for (unsigned int k=0; k<(unsigned int) l.jet_algoPF1_ntk[ijet]; k++) {
      unsigned short trackindex = l.jet_algoPF1_tkind->at(ijet).at(k);
      if (trackindex>=3000) continue; //????
      TLorentzVector TrackP4 = *((TLorentzVector*) l.tk_p4->At(trackindex));
      if (TrackP4.Pt()>1.0) TrackSumP4+=TrackP4;
    }
    if (JetP4.Pt()>MaxJetP4.Pt() && TrackSumP4.Pt()>30) {
      MaxJetP4=JetP4;
      maxjetind=ijet;
      if (TrackSumP4.Pt()>MaxJetTrackP4.Pt()) MaxJetTrackP4=TrackSumP4;
    }
  }
  
  if (maxjetind<0) return false;
  
  if (PADEBUG) cout << "[PADEBUG] : found one good jet" << endl;
  

  vtxAna_.setPairID(maxphoind);

  // now fill the histograms 
  genLevWeight=1.;
  if(cur_type!=0 ) {
    applyGenLevelSmearings(genLevWeight,gP4,l.pu_n,cur_type,genSys,syst_shift);
  }
  evweight  = weight * smeared_pho_weight[maxphoind] * genLevWeight;  
  fillControlPlots(l,photon,MaxJetP4,vtxAna_,maxphoind,evweight);

  return true;

}

// ----------------------------------------------------------------------------------------------------
void PhotonPlusJetVertexAnalysis::fillControlPlots( LoopAll & l, TLorentzVector photon, TLorentzVector jet, HggVertexAnalyzer vtxAna_, int phoindex, float evweight)
{

  // event plots
  l.FillHist("nvtx",0, l.vtx_std_n, evweight);
  
  l.FillHist("photon_pt" , 0 , photon.Pt(),  evweight);
  l.FillHist("photon_eta", 0 , photon.Eta(), evweight);
  l.FillHist("photon_phi", 0 , photon.Phi(), evweight);
  
  l.FillHist("jet_pt", 0  , jet.Pt(),  evweight);
  l.FillHist("jet_eta", 0 , jet.Eta(), evweight);
  l.FillHist("jet_phi", 0 , jet.Phi(), evweight);

  // vertex mva and variables
  // NB vertex 0 is our "true" vertex
  for (int iv = 0 ; iv < l.vtx_std_n; iv++){
    //float dz =  fabs( ( *((TVector3*)l.vtx_std_xyz->At(iv)) -  *((TVector3*)l.vtx_std_xyz->At(0)) ).Z() );  
    //if (dz < 1.){
    if (iv==0){
      // all
      l.FillHist("sumpt2_rv", 0 , vtxAna_.logsumpt2(iv), evweight);
      l.FillHist("ptbal_rv", 0 , vtxAna_.ptbal(iv), evweight);
      l.FillHist("ptasym_rv", 0 , vtxAna_.ptasym(iv), evweight);
      l.FillHist("pulltoconv_rv", 0 , vtxAna_.pulltoconv(iv), evweight);      
      l.FillHist("nconv_rv", 0 , vtxAna_.nconv(iv), evweight);      
      l.FillHist("sumpt2in_rv", 0 , log(vtxAna_.sumpt2in(iv)), evweight);
      l.FillHist("sumpt2out_rv", 0 , log(vtxAna_.sumpt2out(iv)), evweight);
      l.FillHist("nchin_rv", 0 , vtxAna_.nchpho1(iv), evweight);
      l.FillHist("vtxmva_rv", 0 , vtxAna_.mva(iv), evweight);      
      // one good conversion
      if ( vtxAna_.nconv(iv) > 0 ){
	l.FillHist("sumpt2_rv_conv", 0 , vtxAna_.logsumpt2(iv), evweight);
	l.FillHist("ptbal_rv_conv", 0 , vtxAna_.ptbal(iv), evweight);
	l.FillHist("ptasym_rv_conv", 0 , vtxAna_.ptasym(iv), evweight);
	l.FillHist("pulltoconv_rv_conv", 0 , vtxAna_.pulltoconv(iv), evweight);      
	l.FillHist("nconv_rv_conv", 0 , vtxAna_.nconv(iv), evweight);
	l.FillHist("sumpt2in_rv_conv", 0 , log(vtxAna_.sumpt2in(iv)), evweight);
	l.FillHist("sumpt2out_rv_conv", 0 , log(vtxAna_.sumpt2out(iv)), evweight);
	l.FillHist("nchin_rv_conv", 0 , vtxAna_.nchpho1(iv), evweight);
	l.FillHist("vtxmva_rv_conv", 0 , vtxAna_.mva(iv), evweight);      
      }
      // no good conversion, but there are tracks in a cone 0.05 around the photon direction 
      else if (vtxAna_.nchpho1(iv) > 0){
	l.FillHist("sumpt2_rv_tkcone", 0 , vtxAna_.logsumpt2(iv), evweight);
	l.FillHist("ptbal_rv_tkcone", 0 , vtxAna_.ptbal(iv), evweight);
	l.FillHist("ptasym_rv_tkcone", 0 , vtxAna_.ptasym(iv), evweight);
	l.FillHist("pulltoconv_rv_tkcone", 0 , vtxAna_.pulltoconv(iv), evweight);      
	l.FillHist("nconv_rv_tkcone", 0 , vtxAna_.nconv(iv), evweight);
	l.FillHist("sumpt2in_rv_tkcone", 0 , log(vtxAna_.sumpt2in(iv)), evweight);
	l.FillHist("sumpt2out_rv_tkcone", 0 , log(vtxAna_.sumpt2out(iv)), evweight);
	l.FillHist("nchin_rv_tkcone", 0 , vtxAna_.nchpho1(iv), evweight);
	l.FillHist("vtxmva_rv_tkcone", 0 , vtxAna_.mva(iv), evweight);      
      }
    }
    else {
      // all
      l.FillHist("sumpt2_wv", 0 , vtxAna_.logsumpt2(iv), evweight);
      l.FillHist("ptbal_wv", 0 , vtxAna_.ptbal(iv), evweight);
      l.FillHist("ptasym_wv", 0 , vtxAna_.ptasym(iv), evweight);
      l.FillHist("pulltoconv_wv", 0 , vtxAna_.pulltoconv(iv), evweight);      
      l.FillHist("nconv_wv", 0 , vtxAna_.nconv(iv), evweight);      
      l.FillHist("sumpt2in_wv", 0 , log(vtxAna_.sumpt2in(iv)), evweight);
      l.FillHist("sumpt2out_wv", 0 , log(vtxAna_.sumpt2out(iv)), evweight);
      l.FillHist("nchin_wv", 0 , vtxAna_.nchpho1(iv), evweight);
      l.FillHist("vtxmva_wv", 0 , vtxAna_.mva(iv), evweight);      
      // one good conversion
      if ( vtxAna_.nconv(iv) > 0 ){
	l.FillHist("sumpt2_wv_conv", 0 , vtxAna_.logsumpt2(iv), evweight);
	l.FillHist("ptbal_wv_conv", 0 , vtxAna_.ptbal(iv), evweight);
	l.FillHist("ptasym_wv_conv", 0 , vtxAna_.ptasym(iv), evweight);
	l.FillHist("pulltoconv_wv_conv", 0 , vtxAna_.pulltoconv(iv), evweight);      
	l.FillHist("nconv_wv_conv", 0 , vtxAna_.nconv(iv), evweight);
	l.FillHist("sumpt2in_wv_conv", 0 , log(vtxAna_.sumpt2in(iv)), evweight);
	l.FillHist("sumpt2out_wv_conv", 0 , log(vtxAna_.sumpt2out(iv)), evweight);
	l.FillHist("nchin_wv_conv", 0 , vtxAna_.nchpho1(iv), evweight);
	l.FillHist("vtxmva_wv_conv", 0 , vtxAna_.mva(iv), evweight);      
      }
      // no good conversion, but there are tracks in a cone 0.05 around the photon direction 
      else if (vtxAna_.nchpho1(iv) > 0){
	l.FillHist("sumpt2_wv_tkcone", 0 , vtxAna_.logsumpt2(iv), evweight);
	l.FillHist("ptbal_wv_tkcone", 0 , vtxAna_.ptbal(iv), evweight);
	l.FillHist("ptasym_wv_tkcone", 0 , vtxAna_.ptasym(iv), evweight);
	l.FillHist("pulltoconv_wv_tkcone", 0 , vtxAna_.pulltoconv(iv), evweight);      
	l.FillHist("nconv_wv_tkcone", 0 , vtxAna_.nconv(iv), evweight);
	l.FillHist("sumpt2in_wv_tkcone", 0 , log(vtxAna_.sumpt2in(iv)), evweight);
	l.FillHist("sumpt2out_wv_tkcone", 0 , log(vtxAna_.sumpt2out(iv)), evweight);
	l.FillHist("nchin_wv_tkcone", 0 , vtxAna_.nchpho1(iv), evweight);
	l.FillHist("vtxmva_wv_tkcone", 0 , vtxAna_.mva(iv), evweight);      
      }      
    }
  }
  
  
  // event mva 
  int chosenvtx = l.vtx_std_sel;
  float dz =  fabs( ( *((TVector3*)l.vtx_std_xyz->At(chosenvtx)) -  *((TVector3*)l.vtx_std_xyz->At(0)) ).Z() );
  
  l.FillHist("pt", 0 , (photon+jet).Pt(), evweight);

  if (vtxAna_.nconv(chosenvtx) > 0 ) {
    l.FillHist("photon_pt_conv" , 0 , photon.Pt(),  evweight);
    l.FillHist("photon_eta_conv", 0 , photon.Eta(), evweight);
    l.FillHist("pt_conv", 0 , (photon+jet).Pt(), evweight);
  }
  if (vtxAna_.nconv(chosenvtx) == 0 && vtxAna_.nchpho1(chosenvtx) > 0 ) {
    l.FillHist("photon_pt_tkcone" , 0 , photon.Pt(),  evweight);
    l.FillHist("photon_eta_tkcone", 0 , photon.Eta(), evweight);
    l.FillHist("pt_tkcone", 0 , (photon+jet).Pt(), evweight);
  }

  if (dz < 1.) {
    //all
    l.FillHist("evtmva_rv", 0 , l.vtx_std_evt_mva->at(phoindex), evweight);
    l.FillHist("pt_rv", 0 , (photon+jet).Pt(), evweight);
    // one good conversion
    if ( vtxAna_.nconv(chosenvtx) > 0 ){
      l.FillHist("evtmva_rv_conv", 0 , l.vtx_std_evt_mva->at(phoindex), evweight);
      l.FillHist("pt_rv_conv", 0 , (photon+jet).Pt(), evweight);
    }
    // no good conversion but tk in a cone 0.05 around the photon direction
    else if (vtxAna_.nchpho1(chosenvtx) > 0){
      l.FillHist("evtmva_rv_tkcone", 0 , l.vtx_std_evt_mva->at(phoindex), evweight);
      l.FillHist("pt_rv_tkcone", 0 , (photon+jet).Pt(), evweight);
    }      
  }
  else{
    //all
    l.FillHist("evtmva_wv", 0 , l.vtx_std_evt_mva->at(phoindex), evweight);
    // one good conversion
    if ( vtxAna_.nconv(chosenvtx) > 0 ){
      l.FillHist("evtmva_wv_conv", 0 , l.vtx_std_evt_mva->at(phoindex), evweight);
    }
    // no good conversion but tk in a cone 0.05 around the photon direction
    else if (vtxAna_.nchpho1(chosenvtx) > 0){
      l.FillHist("evtmva_wv_tkcone", 0 , l.vtx_std_evt_mva->at(phoindex), evweight);
    }      
  }
  
} 


// ----------------------------------------------------------------------------------------------------
bool PhotonPlusJetVertexAnalysis::SkimEvents(LoopAll& l, int jentry){
    
  return true;

}

// ----------------------------------------------------------------------------------------------------
void PhotonPlusJetVertexAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
}

// ----------------------------------------------------------------------------------------------------
void PhotonPlusJetVertexAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
}

// ----------------------------------------------------------------------------------------------------
bool PhotonPlusJetVertexAnalysis::SelectEventsReduction(LoopAll&, int)
{
    return true;
}
