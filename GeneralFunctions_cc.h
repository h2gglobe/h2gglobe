void LoopAll::GlobeCtIsol(int mode, TLorentzVector* p4, float ptCut, float drCutMin, float drCutMax, Int_t & nIsol, Float_t & ptIsol, Float_t & angle1, Float_t & angle2, Float_t & angle3) {
  nIsol=0;
  ptIsol=0.;
  angle1=10.;
  angle2=10.;
  angle3=10.;

  //must put a track selection

  for (int i=0; i<ct_n; i++) {
    TLorentzVector * tempp4= (TLorentzVector *) ct_p4->At(i);
    if(tempp4->Et()<ptCut) continue;
    double dr=p4->DeltaR(*tempp4);
    if(dr<drCutMin) continue;
    if(dr<angle1) {
      angle3=angle2;
      angle2=angle1;
      angle1=dr;
    }
    else if (dr<angle2) {
      angle3=angle2;
      angle2=dr;
    }
    else if (dr<angle3) {
      angle3=dr;
    }
    if(dr>drCutMax) continue;
    nIsol++;
    ptIsol+=tempp4->Et();
  }
}


int LoopAll::GlobeMatchIsl(TLorentzVector* p4, Float_t & deltaR) {
  deltaR=10.;
  int imatch=-1;

  for (int i=0; i<sc_islbar_n; i++) {
    TLorentzVector * tempp4= (TLorentzVector *) sc_islbar_p4->At(i);
    double dr=p4->DeltaR(*tempp4);
    if(dr<deltaR) {
      deltaR=dr;
      imatch=i;
    }
  }
  if(imatch==-1) {
    cout<<"ERROR GlobeMatchIsl found no match!!!! "<<endl;
  }
  else if(deltaR>0.3) {
    cout<<"Strange, GlobeMatchIsl deltaR="<<deltaR<<" etapho "<<p4->Eta()<<endl;
  }
  return imatch;
}

#include "eIDCuts.h"
std::pair<bool, bool> LoopAll::ElectronId(int index, eIDLevel type) { 

  std::pair<bool, bool> isoIDResult(true, true);

  //TLorentzVector* p4 = (TLorentzVector*)sc_p4->At(el_std_scind[index]);
  TLorentzVector* p4 = (TLorentzVector*)el_std_sc->At(index);
  float eta = fabs(p4->Eta());
  int eb = 0, bin = 0;
  float see = 0;

  if (p4->Et() < 20.) 
    bin = 2;
  else if (p4->Et() > 30.)
    bin = 0;
  else
    bin =1;

  //#ifndef CMSSW3
  //if (eta < 1.479) {
  //see = sc_sieie[el_std_scind[index]];
  //eb = 0;
  //} else {
  //eb = 1; 
  //see = bc_sieie[sc_bcseedind[el_std_scind[index]]];
  //}
  //#else
  if (eta < 1.479) {
    see = el_std_sieiesc[index];
    eb = 0;
  } else {
    eb = 1; 
    see = el_std_sieie[index];
  }
  //#endif

  float eseedopincor = el_std_eseedopin[index] + el_std_fbrem[index];

  if(el_std_fbrem[index]<0) 
    eseedopincor = el_std_eseedopin[index]; 

  //#ifndef CMSSW3
  //float sip = sipCalculator(index);
  //#else
  float sip = fabs(el_std_ip_gsf[index]);
  //#endif

  int cat = ElectronClassification(index);
  
  float corr_tk_iso   = el_std_tkiso03[  index];
  float corr_ecal_iso = el_std_ecaliso04[index];
  float corr_hcal_iso = el_std_hcaliso04[index];  

  corr_tk_iso   = corr_tk_iso  *pow(40/p4->Et(), 2); 
  corr_ecal_iso = corr_ecal_iso*pow(40/p4->Et(), 2);
  corr_hcal_iso = corr_hcal_iso*pow(40/p4->Et(), 2);
  
  if ((corr_tk_iso > cutisotk[bin][type][cat]) ||
      (corr_ecal_iso > cutisoecal[bin][type][cat]) ||
      (corr_hcal_iso > cutisohcal[bin][type][cat]))

    isoIDResult.first = false;
  
  if (el_std_fbrem[index] < -2) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (el_std_hoe[index]  > cuthoe[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (see > cutsee[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (fabs(el_std_dphiin[index]) > cutdphi[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;  
  }                            

  if (fabs(el_std_detain[index]) > cutdeta[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (eseedopincor < cuteopin[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }                      

  if (sip > cutip[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (el_std_hp_expin[index]  > cutmishits[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  return isoIDResult;
}

int LoopAll::ElectronClassification(int index) {
  TLorentzVector* p4 = (TLorentzVector*) el_std_sc->At(index);

  int cat = -1;
  float eta = fabs(p4->Eta());

  if (eta < 1.479) {       // BARREL
    if(el_std_fbrem[index]<0.12)
      cat=1;
    else if (el_std_eopin[index] < 1.2 && el_std_eopin[index] > 0.9) 
      cat=0;
    else 
      cat=2;
  } else {                // ENDCAP
    if(el_std_fbrem[index]<0.2)
      cat=4;
    else if (el_std_eopin[index] < 1.22 && el_std_eopin[index] > 0.82) 
      cat=3;
    else 
      cat=5;
  }

  return cat;
}

Float_t LoopAll::sipCalculator(int index) {

  Float_t ip = 0;

  if (el_std_tkind[index] != -1) {
    TLorentzVector* tk = (TLorentzVector*)tk_p4->At(el_std_tkind[index]);
    TVector3* my_tk_pos = (TVector3*)tk_vtx_pos->At(el_std_tkind[index]);

    // FIXME to handle the case of multiple vertices

    if (vtx_std_n != 0) {
      TVector3* my_vtx_pos = (TVector3*)vtx_std_xyz->At(0);

      // this is d0 "corrected" for the vertex...
      ip = fabs((-(my_tk_pos->X()-my_vtx_pos->X())*tk->Y()+(my_tk_pos->Y()-my_vtx_pos->Y()) * tk->X())/tk->Pt());
    } else {
      ip = fabs((-(my_tk_pos->X())*tk->Y()+(my_tk_pos->Y()) * tk->X())/tk->Pt());
    }

  }

  return ip;
}

void LoopAll::eIDInfo(Int_t index, Int_t& iso_result, Int_t& id_result, Int_t eIDMaxLevel) {

  iso_result = 0;
  id_result = 0;

  // FIXME add GetEntry functions
  for(Int_t i=0; i<eIDMaxLevel; ++i) {
    std::pair<bool, bool> result = ElectronId(index, LoopAll::eIDLevel(i));

    if (result.first) 
      iso_result = i;

    if (result.second)
      id_result = i;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
// Vertex Analysis
// ---------------------------------------------------------------------------------------------------------------------------------------------
class GlobeVertexInfo : public VertexInfoAdapter
{
public:
	GlobeVertexInfo(LoopAll &);
	
	virtual int nvtx() const    { return lo_.vtx_std_n; };
	virtual int ntracks() const { return lo_.tk_n; };
	
	virtual float tkpx(int ii) const { return ((TVector3*)lo_.tk_p4->At(ii))->Px(); };
	virtual float tkpy(int ii) const { return ((TVector3*)lo_.tk_p4->At(ii))->Py(); };
	virtual float tkpz(int ii) const { return ((TVector3*)lo_.tk_p4->At(ii))->Pz(); };
	
	virtual float tkPtErr(int ii) const { return lo_.tk_pterr[ii]; };
	virtual int   tkVtxId(int ii) const { return -1; }; // FIXME

	virtual float tkWeight(int ii) const { return -1.;}  // FIXME
	
	virtual float vtxx(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->X(); };
	virtual float vtxy(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Y(); };
	virtual float vtxz(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Z(); };

	virtual float tkd0(int ii) const { return 0.; } // FIXME
	virtual float tkd0Err(int ii) const { return 0.; };  // FIXME

	virtual float tkdz(int ii) const { return 0.; };  // FIXME
	virtual float tkdzErr(int ii) const { return 0.; };  // FIXME

	virtual bool tkIsHighPurity(int ii) const { return ( lo_.tk_quality[ii] & (1<<2) ) >> 2; };

	virtual ~GlobeVertexInfo();
	
private:
	LoopAll & lo_;
};


// ---------------------------------------------------------------------------------------------------------------------------------------------
GlobeVertexInfo::GlobeVertexInfo(LoopAll & lo) : lo_(lo) {};

// ---------------------------------------------------------------------------------------------------------------------------------------------
GlobeVertexInfo::~GlobeVertexInfo() {};


// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::vertexAnalysis(int p1, int p2)
{
	GlobeVertexInfo vinfo(*this);
	
	PhotonInfo pho1(*((TVector3*)pho_calopos->At(p1)),((TLorentzVector*)pho_p4->At(p1))->Energy()),pho2(*((TVector3*)pho_calopos->At(p2)),((TLorentzVector*)pho_p4->At(p2))->Energy());

	HggVertexAnalyzer vAna(vtxAlgoParams,vtx_std_n);
	
	vAna.analyze(vinfo,pho1,pho2);
	
	std::vector<int> rankprod = vAna.rankprod(vtxVarNames);
	cout << "\n\nRanks product" << endl;
	cout << "best vertex " << rankprod[0] << endl;
	for(int ii=0; ii<vtx_std_n; ++ii) {
		int vtxrank = find(rankprod.begin(), rankprod.end(), ii) - rankprod.begin();
		cout << "vertx " << ii << " rank " << vtxrank << " " << vAna.ptbal(ii) << " " << vAna.ptasym(ii) << " " << vAna.logsumpt2(ii) << endl;
		
	}
	
}


