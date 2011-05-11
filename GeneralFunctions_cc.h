// ---------------------------------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------------------------------------------------------------------------
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
	
	virtual bool hasVtxTracks() const { return true; }
	virtual const unsigned short * vtxTracks(int ii) const { return &(*lo_.vtx_std_tkind)[ii][0]; };
	virtual int vtxNTracks(int ii) const { return lo_.vtx_std_ntks[ii]; };
	virtual const float * vtxTkWeights(int ii) const { return &(*lo_.vtx_std_tkweight)[ii][0]; };

	virtual float tkpx(int ii) const { return ((TLorentzVector*)lo_.tk_p4->At(ii))->Px(); };
	virtual float tkpy(int ii) const { return ((TLorentzVector*)lo_.tk_p4->At(ii))->Py(); };
	virtual float tkpz(int ii) const { return ((TLorentzVector*)lo_.tk_p4->At(ii))->Pz(); };
	
	virtual float tkPtErr(int ii) const { return lo_.tk_pterr[ii]; };
	virtual int   tkVtxId(int ii) const { return -1; };

	virtual float tkWeight(int ii, int jj) const { return (*lo_.vtx_std_tkweight)[jj][ii]; };
	
	virtual float vtxx(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->X(); };
	virtual float vtxy(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Y(); };
	virtual float vtxz(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Z(); };

	virtual float tkd0(int ii, int jj) const { return 0.; } // FIXME
	virtual float tkd0Err(int ii, int jj) const { return 0.; };  // FIXME

	virtual float tkdz(int ii, int jj) const { return 0.; };  // FIXME
	virtual float tkdzErr(int ii, int jj) const { return 0.; };  // FIXME

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
void LoopAll::vertexAnalysis(HggVertexAnalyzer & vtxAna, int p1, int p2)
{
        GlobeVertexInfo vinfo(*this); 
	PhotonInfo
	  pho1(p1,*((TVector3*)pho_calopos->At(p1)),((TLorentzVector*)pho_p4->At(p1))->Energy()),
	  pho2(p2,*((TVector3*)pho_calopos->At(p2)),((TLorentzVector*)pho_p4->At(p2))->Energy());
	vtxAna.analyze(vinfo,pho1,pho2);
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> LoopAll::vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, int p1, int p2, std::vector<std::string> & vtxVarNames)
{
	assert( p1 == vtxAna.pho1() && p2 == vtxAna.pho2() );

	// preselect vertices : all vertices
        std::vector<int> preselAll;
        for(int i=0; i<vtx_std_n ; i++) {
          preselAll.push_back(i); 
        }

	// conversions infos
	PhotonInfo pho1(p1,
			*((TVector3*)pho_calopos->At(p1)),
			*((TVector3*) bs_xyz),
			*((TVector3*) pho_conv_vtx->At(p1)),
			((TLorentzVector*)pho_p4->At(p1))->Energy(),
			pho_isEB[p1],
			pho_conv_ntracks[p1],
			pho_conv_validvtx[p1],
			pho_conv_chi2_probability[p1] ,
			pho_conv_eoverp[p1]
			);
	

	PhotonInfo pho2(p2,
			*((TVector3*)pho_calopos->At(p2)),
			*((TVector3*) bs_xyz),
			*((TVector3*) pho_conv_vtx->At(p2)),
			((TLorentzVector*)pho_p4->At(p2))->Energy(),
			pho_isEB[p2],
			pho_conv_ntracks[p2],
			pho_conv_validvtx[p2],
			pho_conv_chi2_probability[p2] ,
			pho_conv_eoverp[p2]
			);
		
        float zconv = 0; 
        float dzconv = 0;
        std::vector<int> preselConv;


        if ( (pho_r9[p1] <0.93 || pho_r9[p2] <0.93) && (pho1.isAConversion() || pho2.isAConversion()) )  {
	  
          if (pho1.isAConversion()  && !pho2.isAConversion() ){
            zconv  = vtxAnaFromConv.vtxZ(pho1);
            dzconv = vtxAnaFromConv.vtxdZ(pho1);
          }
	  
          if (pho2.isAConversion() && !pho1.isAConversion()){
            zconv  = vtxAnaFromConv.vtxZ(pho2);
            dzconv = vtxAnaFromConv.vtxdZ(pho2);
          }
	  
          if ( pho1.isAConversion() && pho2.isAConversion()){
            float z1  = vtxAnaFromConv.vtxZ(pho1);
            float dz1 = vtxAnaFromConv.vtxdZ(pho1);
            
            float z2  = vtxAnaFromConv.vtxZ(pho2);
            float dz2 = vtxAnaFromConv.vtxdZ(pho2);
            
            zconv  = sqrt ( 1./(1./dz1/dz1 + 1./dz2/dz2 )*(z1/dz1/dz1 + z2/dz2/dz2) ) ;  // weighted average
            dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
          }
	  
	  // preselect vertices : only vertices in a window zconv +/- dzconv

	  for(int i=0; i < vtx_std_n; i++) {
	    TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(i);
	    if ( fabs(zconv - vtxpos->Z() ) < dzconv ) 
              preselConv.push_back(i); 
          }
	  
        }

	// preselection 
	if ( preselConv.size()==0 ) 
          vtxAna.preselection(preselAll);
        else 
          vtxAna.preselection(preselConv);

	std::vector<int> rankprod = vtxAna.rankprod(vtxVarNames);
	//// cout << "\n\nRanks product" << endl;
	//// cout << "best vertex " << rankprod[0] << endl;
	//// for(int ii=0; ii<vtx_std_n; ++ii) {
	//// 	int vtxrank = find(rankprod.begin(), rankprod.end(), ii) - rankprod.begin();
	//// 	cout << "vertx " << ii << " rank " << vtxrank << " " << vtxAna.ptbal(ii) << " " << vtxAna.ptasym(ii) << " " << vtxAna.logsumpt2(ii) << endl;
	//// }
	
	return rankprod;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector LoopAll::get_pho_p4(int ipho, int ivtx)
{
	PhotonInfo p(ipho, *((TVector3*)pho_calopos->At(ipho)),((TLorentzVector*)pho_p4->At(ipho))->Energy());
	TVector3 * vtx = (TVector3*) vtx_std_xyz->At(ivtx);
	return p.p4( vtx->X(), vtx->Y(), vtx->Z() );
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::set_pho_p4(int ipho, int ivtx)
{
	PhotonInfo p(ipho, *((TVector3*)pho_calopos->At(ipho)),((TLorentzVector*)pho_p4->At(ipho))->Energy());
	TVector3 * vtx = (TVector3*) vtx_std_xyz->At(ivtx);
	TLorentzVector p4 = p.p4( vtx->X(), vtx->Y(), vtx->Z() );
	*((TLorentzVector*)pho_p4->At(ipho)) = p4;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::FillCICInputs()
{
	pho_tkiso_recvtx_030_002_0000_10_01->clear(); pho_tkiso_recvtx_030_002_0000_10_01->resize(pho_n,std::vector<float>(vtx_std_n,0.));
	
	for(int ipho=0;ipho<pho_n;++ipho){
		TLorentzVector * phop4 = (TLorentzVector*)pho_p4->At(ipho);
		pho_tkiso_badvtx_040_002_0000_10_01[ipho] = WorstSumTrackPtInCone(phop4, 0,0, 0.40, 0.02, 0.0, 1.0, 0.1);
		pho_drtotk_25_99[ipho] = DeltaRToTrack(ipho, vtx_std_sel, 2.5, 99.);
		
		for(int ivtx=0;ivtx<vtx_std_n;++ivtx) {
			TLorentzVector p4 = get_pho_p4( ipho, ivtx );
			(*pho_tkiso_recvtx_030_002_0000_10_01)[ipho][ivtx] = SumTrackPtInCone(phop4, ivtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::FillCIC()
{
	pho_passcuts_lead->clear(); pho_passcuts_lead->resize( pho_n, std::vector<UInt_t>(phoNCUTLEVELS,0) ); 
	pho_cutlevel_lead->clear(); pho_cutlevel_lead->resize( pho_n, 0 );
	pho_passcuts_sublead->clear(); pho_passcuts_sublead->resize( pho_n, std::vector<UInt_t>(phoNCUTLEVELS,0) ); 
	pho_cutlevel_sublead->clear(); pho_cutlevel_sublead->resize( pho_n, 0 );
	
	std::vector<std::vector<bool> > passcut_lead, passcut_sublead;
	for(int ipho=0;ipho<pho_n;++ipho){
		int level_lead = PhotonCiCSelectionLevel(ipho, passcut_lead, 0);
		int level_sublead = PhotonCiCSelectionLevel(ipho, passcut_sublead, 1);
		(*pho_cutlevel_lead)[ipho] = level_lead;
		(*pho_cutlevel_sublead)[ipho] = level_sublead;
		for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {
			UInt_t leadw=0, subleadw=0;
			for(size_t icut=0; icut<passcut_lead.size(); ++icut) {
				leadw |= ( (!passcut_lead[iCUTLEVEL][icut] & 0x1) << icut);
			}
			(*pho_passcuts_lead)[ipho][iCUTLEVEL] = leadw;
			for(size_t icut=0; icut<passcut_sublead.size(); ++icut) {
				subleadw |= ( (!passcut_sublead[iCUTLEVEL][icut] & 0x1) << icut);
			}
			(*pho_passcuts_sublead)[ipho][iCUTLEVEL] = subleadw;
			
		}
	}
}


// CiC SELECTION CODE BEGIN - SSIMON
// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * allcuts_lead, float * allcuts_sublead) {

  //thresholds are in this order below
  // isosumoet[]
  // isosumoetbad[]
  // trkisooetom[]
  // sieie[]
  // hovere[]
  // r9[]
  // drtotk_25_99[]
  // pixel[]

// phoLOOSE       - sob value=0.000307692    - iteration 8 - eff=0.95544 fake=0.138351
// phoMEDIUM      - sob value=0.000615385    - iteration 6 - eff=0.93371 fake=0.0992173
// phoTIGHT       - sob value=0.00123077     - iteration 6 - eff=0.901978 fake=0.0720423
// phoSUPERTIGHT  - sob value=0.00246154     - iteration 6 - eff=0.844751 fake=0.0426593
// phoHYPERTIGHT1 - sob value=0.00492308     - iteration 6 - eff=0.781929 fake=0.0277505
// phoHYPERTIGHT2 - sob value=0.00961538     - iteration 6 - eff=0.692504 fake=0.0170745
// phoHYPERTIGHT3 - sob value=0.0192308      - iteration 6 - eff=0.560266 fake=0.00916371
// phoHYPERTIGHT4 - sob value=0.0384615      - iteration 6 - eff=0.409147 fake=0.00415847

  const unsigned int ncuts = 8;
  const unsigned int ncategories = 6;
  switch(cutlevel) {
    case(phoLOOSE) : {
                    float allcuts_temp_lead[] = { 
                      13.2932,     11.6658,     10.3763,     10.8621,     10.2079,     8.73226,
                      31.5827,     21.7584,     14.4519,      18.342,     14.2522,     12.2147,
                      7.27839,     5.40834,      3.2986,     5.35507,     5.07699,     3.12837,
                      0.0108999,   0.0108493,    0.010071,   0.0295929,   0.0288815,   0.0275799,
                      0.0925256,   0.0972961,   0.0953947,   0.0908323,   0.0812638,   0.0485916,
                      0.95000,     0.89999,    0.260709,    0.9500,    0.90000,    0.395383,
                      98.0036,     98.0038,   0.0132181,  0.00642242,  0.00687745,   0.0121392,
                      1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                    float allcuts_temp_sublead[] = { 
                      13.2932,     11.6658,     10.3763,     10.8621,     10.2079,     8.73226,
                      31.5827,     21.7584,     14.4519,      18.342,     14.2522,     12.2147,
                      7.27839,     5.40834,      3.2986,     5.35507,     5.07699,     3.12837,
                      0.0108999,   0.0108493,    0.010071,   0.0295929,   0.0288815,   0.0275799,
                      0.0925256,   0.0972961,   0.0953947,   0.0908323,   0.0812638,   0.0485916,
                      0.95000,     0.89999,    0.260709,    0.9500,    0.90000,    0.395383,
                      98.0036,     98.0038,   0.0132181,  0.00642242,  0.00687745,   0.0121392,
                      1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                    for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                      allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                  } break;
    case(phoMEDIUM) : {
                     float allcuts_temp_lead[] = {  
                       11.9631,     10.4161,     9.29982,     10.8621,       9.734,      8.4268,
                       20.8853,     16.4926,     13.5553,     15.2879,     11.4552,     10.7918,
                       6.32821,     5.39826,     3.06333,     3.94555,     3.65483,     2.31875,
                       0.0108495,    0.010793,  0.00976376,   0.0295929,   0.0287852,    0.027411,
                       0.0924405,   0.0972961,   0.0644553,   0.0908323,   0.0765037,   0.0485916,
                       0.95000,    0.9000,    0.262508,    0.9500,    0.90005,    0.444622,
                       98.9721,     98.9721,   0.0164165,  0.00642242,  0.00940864,   0.0126402,
                       1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                     float allcuts_temp_sublead[] = {  
                       11.9631,     10.4161,     9.29982,     10.8621,       9.734,      8.4268,
                       20.8853,     16.4926,     13.5553,     15.2879,     11.4552,     10.7918,
                       6.32821,     5.39826,     3.06333,     3.94555,     3.65483,     2.31875,
                       0.0108495,    0.010793,  0.00976376,   0.0295929,   0.0287852,    0.027411,
                       0.0924405,   0.0972961,   0.0644553,   0.0908323,   0.0765037,   0.0485916,
                       0.95000,    0.9000,    0.262508,    0.9500,    0.90005,    0.444622,
                       98.9721,     98.9721,   0.0164165,  0.00642242,  0.00940864,   0.0126402,
                       1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                     for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                       allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                   } break;
    case(phoTIGHT) : {
                    float allcuts_temp_lead[] = { 
                      11.0997,     9.53812,     8.73053,     10.5643,     8.20721,      7.5888,
                      18.7516,     13.8868,     11.7356,     13.3007,     10.1418,      9.9876,
                      5.83477,     4.64656,     2.55216,     3.45855,     3.54738,     2.04094,
                      0.0108495,   0.0107461,  0.00960384,   0.0295009,   0.0279728,   0.0269579,
                      0.0924405,   0.0970314,    0.023329,   0.0688847,   0.0763571,   0.0316924,
                      0.95000,    0.9000,    0.283541,    0.9500,    0.9000,     0.45754,
                      98.9992,     98.9992,   0.0200119,  0.00642242,  0.00940864,   0.0144829,
                      1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                    float allcuts_temp_sublead[] = { 
                      11.0997,     9.53812,     8.73053,     10.5643,     8.20721,      7.5888,
                      18.7516,     13.8868,     11.7356,     13.3007,     10.1418,      9.9876,
                      5.83477,     4.64656,     2.55216,     3.45855,     3.54738,     2.04094,
                      0.0108495,   0.0107461,  0.00960384,   0.0295009,   0.0279728,   0.0269579,
                      0.0924405,   0.0970314,    0.023329,   0.0688847,   0.0763571,   0.0316924,
                      0.95000,    0.9000,    0.283541,    0.9500,    0.9000,     0.45754,
                      98.9992,     98.9992,   0.0200119,  0.00642242,  0.00940864,   0.0144829,
                      1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                    for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                      allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                  } break;
    case(phoSUPERTIGHT) : {
                         float allcuts_temp_lead[] = { 
                           10.0747,     8.81289,     8.19314,     8.14119,     7.25452,     7.32875,
                           14.8691,     13.0847,     10.1978,     11.9686,     9.29623,     9.55426,
                           3.65654,     3.07026,     2.33651,     2.79337,     3.16278,      1.4586,
                           0.0102523,   0.0103579,  0.00952823,   0.0295009,    0.026189,    0.025918,
                           0.0924405,   0.0970314,    0.023329,   0.0688847,   0.0647023,    0.023614,
                           0.95000,    0.9000,    0.298177,    0.9500,    0.9000,    0.459521,
                           98,          98,   0.0214232,  0.00642242,     13.6939,   0.0161272,
                           1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                         float allcuts_temp_sublead[] = { 
                           10.0747,     8.81289,     8.19314,     8.14119,     7.25452,     7.32875,
                           14.8691,     13.0847,     10.1978,     11.9686,     9.29623,     9.55426,
                           3.65654,     3.07026,     2.33651,     2.79337,     3.16278,      1.4586,
                           0.0102523,   0.0103579,  0.00952823,   0.0295009,    0.026189,    0.025918,
                           0.0924405,   0.0970314,    0.023329,   0.0688847,   0.0647023,    0.023614,
                           0.95000,    0.9000,    0.298177,    0.9500,    0.9000,    0.459521,
                           98,          98,   0.0214232,  0.00642242,     13.6939,   0.0161272,
                           1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                         for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                           allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                       } break;
    case(phoHYPERTIGHT1) : {
                          float allcuts_temp_lead[] = { 
                            9.76772,     7.96695,     7.20682,     8.01842,     6.73697,     6.67724,
                            12.8069,     11.8664,     9.45415,     11.4739,     8.05775,      9.1949,
                            3.45884,     3.07026,      1.4997,     1.71704,     3.16278,     1.24066,
                            0.0101858,   0.0102206,  0.00931955,   0.0276087,   0.0256188,   0.0248065,
                            0.0910159,   0.0970314,    0.023329,   0.0578315,   0.0584827,    0.023614,
                            0.95000,    0.9000,    0.359721,    0.9500,    0.9000,    0.459521,
                            98,          98,     1.88292,     5.14027,     92.3666,     1.85543,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          float allcuts_temp_sublead[] = { 
                            9.76772,     7.96695,     7.20682,     8.01842,     6.73697,     6.67724,
                            12.8069,     11.8664,     9.45415,     11.4739,     8.05775,      9.1949,
                            3.45884,     3.07026,      1.4997,     1.71704,     3.16278,     1.24066,
                            0.0101858,   0.0102206,  0.00931955,   0.0276087,   0.0256188,   0.0248065,
                            0.0910159,   0.0970314,    0.023329,   0.0578315,   0.0584827,    0.023614,
                            0.95000,    0.9000,    0.359721,    0.9500,    0.9000,    0.459521,
                            98,          98,     1.88292,     5.14027,     92.3666,     1.85543,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                            allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                        } break;
    case(phoHYPERTIGHT2) : {
                          float allcuts_temp_lead[] = { 
                            8.42342,     7.32305,     6.70688,     6.87994,     6.73697,     6.28105,
                            12.6492,     9.62967,     8.60529,     8.53858,     6.74444,     8.31702,
                            2.93957,     2.67696,      1.2676,     1.71704,     2.67272,     1.06582,
                            0.0100159,   0.0102206,  0.00931955,   0.0276087,   0.0254186,   0.0231607,
                            0.0571565,   0.0918075,   0.0193085,   0.0565089,   0.0584827,    0.023614,
                            0.95000,    0.900,    0.369148,    0.950,    0.90,    0.459521,
                            98,          98,     1.96832,     96.3725,     92.3666,     1.85543,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          float allcuts_temp_sublead[] = { 
                            8.42342,     7.32305,     6.70688,     6.87994,     6.73697,     6.28105,
                            12.6492,     9.62967,     8.60529,     8.53858,     6.74444,     8.31702,
                            2.93957,     2.67696,      1.2676,     1.71704,     2.67272,     1.06582,
                            0.0100159,   0.0102206,  0.00931955,   0.0276087,   0.0254186,   0.0231607,
                            0.0571565,   0.0918075,   0.0193085,   0.0565089,   0.0584827,    0.023614,
                            0.95000,    0.900,    0.369148,    0.950,    0.90,    0.459521,
                            98,          98,     1.96832,     96.3725,     92.3666,     1.85543,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                            allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                        } break;
    case(phoHYPERTIGHT3) : {
                          float allcuts_temp_lead[] = { 
                            7.54771,     6.87163,     6.13022,     6.87994,     5.81362,     6.14728,
                            11.3148,      8.8786,     7.98267,     6.85881,     5.40019,     8.01959,
                            2.3804,     1.64676,     1.14244,     1.71704,     1.81568,     1.05839,
                            0.009886,  0.00979577,  0.00900805,   0.0276087,   0.0254186,   0.0218532,
                            0.0492746,   0.0233285,   0.0193085,   0.0405891,   0.0584827,    0.023614,
                            0.95000,    0.90,    0.369148,    0.950,    0.90,    0.459521,
                            98,          98,     5.03113,     97.4235,     92.3666,     8.32375,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          float allcuts_temp_sublead[] = { 
                            7.54771,     6.87163,     6.13022,     6.87994,     5.81362,     6.14728,
                            11.3148,      8.8786,     7.98267,     6.85881,     5.40019,     8.01959,
                            2.3804,     1.64676,     1.14244,     1.71704,     1.81568,     1.05839,
                            0.009886,  0.00979577,  0.00900805,   0.0276087,   0.0254186,   0.0218532,
                            0.0492746,   0.0233285,   0.0193085,   0.0405891,   0.0584827,    0.023614,
                            0.95000,    0.90,    0.369148,    0.950,    0.90,    0.459521,
                            98,          98,     5.03113,     97.4235,     92.3666,     8.32375,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                            allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                        } break;
    case(phoHYPERTIGHT4) : {
                          float allcuts_temp_lead[] = { 
                            6.47675,     6.87163,     5.69614,     5.91569,     5.81362,     6.14728,
                            9.48083,     8.41617,     7.47837,     6.77279,     5.40019,     6.98873,
                            1.19423,     1.45099,    0.791635,     1.71704,     1.08941,    0.635031,
                            0.009886,   0.0095335,  0.00867298,   0.0261185,   0.0240073,   0.0215452,
                            0.0492746,   0.0164099,  0.00397353,   0.0346242,   0.0350896,   0.0184793,
                            0.9500,    0.90,    0.480906,    0.95,     0.90,    0.606315,
                            98,          98,     96.3695,     97.4235,      96.612,     8.32375,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          float allcuts_temp_sublead[] = { 
                            6.47675,     6.87163,     5.69614,     5.91569,     5.81362,     6.14728,
                            9.48083,     8.41617,     7.47837,     6.77279,     5.40019,     6.98873,
                            1.19423,     1.45099,    0.791635,     1.71704,     1.08941,    0.635031,
                            0.009886,   0.0095335,  0.00867298,   0.0261185,   0.0240073,   0.0215452,
                            0.0492746,   0.0164099,  0.00397353,   0.0346242,   0.0350896,   0.0184793,
                            0.9500,    0.90,    0.480906,    0.95,     0.90,    0.606315,
                            98,          98,     96.3695,     97.4235,      96.612,     8.32375,
                            1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                          for(int i=0;i!=ncuts*ncategories;++i) { allcuts_lead[i]=allcuts_temp_lead[i];
                            allcuts_sublead[i]=allcuts_temp_sublead[i]; }
                        } break;
    default:std::cout << "UNKNOWN phoCiCIDLevel: " << cutlevel << std::endl;

  }
}



// ---------------------------------------------------------------------------------------------------------------------------------------------
std::pair<int,int> LoopAll::DiphotonCiCSelection( phoCiCIDLevel LEADCUTLEVEL, phoCiCIDLevel SUBLEADCUTLEVEL, Float_t leadPtMin, Float_t subleadPtMin, bool applyPtoverM) {

  //rho=0;// CAUTION SETTING RHO TO 0 FOR 2010 DATA FILES (RHO ISN'T IN THESE FILES)
  int selected_lead_index = -1;
  int selected_sublead_index = -1;
  float selected_lead_pt = -1;
  float selected_sublead_pt = -1;

  std::vector<std::vector<bool> > ph_passcut;
  for(int ipho=0;ipho!=pho_n;++ipho) {
    TLorentzVector * iphop4 = (TLorentzVector*)pho_p4->At(ipho);
    if(iphop4->Et() < leadPtMin || fabs(iphop4->Eta()) > 2.5)continue;

    if(PhotonCiCSelectionLevel(ipho, ph_passcut, 0) < LEADCUTLEVEL)continue;

    for(int iipho=0;iipho!=pho_n;++iipho) {
      if(iipho == ipho)continue;
      TLorentzVector * iiphop4 = (TLorentzVector*)pho_p4->At(iipho);
      if(iiphop4->Et() < subleadPtMin || fabs(iiphop4->Eta()) > 2.5)continue;
      if(iiphop4->Et() > iphop4->Et())continue;
      float m_gamgam = (*iphop4+*iiphop4).M();
      float L_ptom = iphop4->Et()/m_gamgam;
      float S_ptom = iiphop4->Et()/m_gamgam;
      if(applyPtoverM && (L_ptom < 0.33 || S_ptom<0.25))continue;

      if(PhotonCiCSelectionLevel(iipho, ph_passcut, 1) < SUBLEADCUTLEVEL )continue;
      // if here, diphoton passed all cuts.
      //std::cout << "FOUND DIPHOTON" << std::endl;
      if( (iphop4->Et()>selected_lead_pt) || (ipho==selected_lead_index&&iiphop4->Et()>selected_sublead_pt) ) {
	      selected_lead_pt = iphop4->Et();
	      selected_sublead_pt = iiphop4->Et();
	      selected_lead_index = ipho;
	      selected_sublead_index = iipho;
      }
      
    }// end photon loop (iipho), aka sublead
  }// end photon loop (ipho), aka lead

  std::pair<int,int> dipho_inds(selected_lead_index,selected_sublead_index);
  return dipho_inds;

}


// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::PhotonCiCSelectionLevel( int photon_index, std::vector<std::vector<bool> > & ph_passcut, int doSublead) {

  int cutlevelpassed = -1;

  int n_r9_categories = 3;
  int n_eta_categories = 2;
  int photon_category = PhotonCategory(photon_index,n_r9_categories,n_eta_categories);

  TLorentzVector * phop4 = (TLorentzVector*)pho_p4->At(photon_index);

  float val_tkiso = (*pho_tkiso_recvtx_030_002_0000_10_01)[photon_index][vtx_std_sel];
  float val_ecaliso = pho_ecalsumetconedr03[photon_index];
  float val_hcaliso = pho_hcalsumetconedr04[photon_index];
  float val_ecalisobad = pho_ecalsumetconedr04[photon_index];
  float val_hcalisobad = pho_hcalsumetconedr04[photon_index];
  float val_tkisobad = pho_tkiso_badvtx_040_002_0000_10_01[photon_index];
  float val_sieie = pho_sieie[photon_index];
  float val_hoe = pho_hoe[photon_index];
  float val_r9 = pho_r9[photon_index];
  float val_drtotk_25_99 = pho_drtotk_25_99[photon_index];
  float val_pixel = (float)pho_haspixseed[photon_index];

  float rhofacbad=0.40, rhofac=0.05;
  float val_isosumoet=(val_tkiso+val_ecaliso+val_hcaliso+5-rho*rhofac)*50./phop4->Et();
  float val_isosumoetbad=(val_tkisobad+val_ecalisobad+val_hcalisobad+7-rho*rhofacbad)*50./phop4->Et();
  float val_trkisooet=(val_tkiso)*50./phop4->Et();

  ph_passcut.clear();
  ph_passcut.resize(phoNCUTLEVELS,std::vector<bool>(8,true) );
  if(!doSublead) {
    for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {
      ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <   cut_lead_isosumoet[iCUTLEVEL][photon_category]     );
      ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <   cut_lead_isosumoetbad[iCUTLEVEL][photon_category]  );
      ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <   cut_lead_trkisooet[iCUTLEVEL][photon_category]     );
      ph_passcut[iCUTLEVEL][3] = (val_sieie            <   cut_lead_sieie[iCUTLEVEL][photon_category]         );
      ph_passcut[iCUTLEVEL][4] = (val_hoe              <   cut_lead_hovere[iCUTLEVEL][photon_category]        );
      ph_passcut[iCUTLEVEL][5] = (val_r9             >     cut_lead_r9[iCUTLEVEL][photon_category]            );// gt cut
      ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >     cut_lead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
      ph_passcut[iCUTLEVEL][7] = (val_pixel            <   cut_lead_pixel[iCUTLEVEL][photon_category]         );
      bool ph_passcut_all = true;
      for(int icut=0;icut!=8;++icut) {
	ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
      }
      if(ph_passcut_all) {
	if( cutlevelpassed != iCUTLEVEL - 1 ) {
	  std::cerr << "photon " << photon_index << " (category " << photon_category << ") in run/event " << run << "/" << event << " passed CiC cut level " 
		    << iCUTLEVEL << " but not "  << iCUTLEVEL - 1 << ". Did you load your cut values correctly? "<< std::endl;
	  /// assert( 0 );
	}
	cutlevelpassed=iCUTLEVEL;
      }
    }
  } else if(doSublead) {
    for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {
      ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <   cut_sublead_isosumoet[iCUTLEVEL][photon_category]     );
      ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <   cut_sublead_isosumoetbad[iCUTLEVEL][photon_category]  );
      ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <   cut_sublead_trkisooet[iCUTLEVEL][photon_category]     );
      ph_passcut[iCUTLEVEL][3] = (val_sieie            <   cut_sublead_sieie[iCUTLEVEL][photon_category]         );
      ph_passcut[iCUTLEVEL][4] = (val_hoe              <   cut_sublead_hovere[iCUTLEVEL][photon_category]        );
      ph_passcut[iCUTLEVEL][5] = (val_r9             >     cut_sublead_r9[iCUTLEVEL][photon_category]            );// gt cut
      ph_passcut[iCUTLEVEL][6] = (val_drtotk_25_99   >     cut_sublead_drtotk_25_99[iCUTLEVEL][photon_category]  );// gt cut
      ph_passcut[iCUTLEVEL][7] = (val_pixel            <   cut_sublead_pixel[iCUTLEVEL][photon_category]         );
      bool ph_passcut_all = true;
      for(int icut=0;icut!=8;++icut) {
	ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
      }
      if(ph_passcut_all) {
	if( cutlevelpassed != iCUTLEVEL - 1 ) {
	  std::cerr << "photon " << photon_index << " (category " << photon_category << ") in run/event " << run << "/" << event << " passed CiC cut level " 
		    << iCUTLEVEL << " but not "  << iCUTLEVEL - 1 << ". Did you load your cut values correctly? " << std::endl;
	    assert( 0 );
	}
	cutlevelpassed=iCUTLEVEL;
      }
    }
  }

  return cutlevelpassed;

}


// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::DeltaRToTrack(Int_t photonind, Int_t vtxind, Float_t PtMin, Float_t dzmax, Float_t dxymax, int maxlosthits){
  if(LDEBUG)std::cout << "DeltaRToTrack BEGIN" << std::endl;
  int elind = -1;
  float eldr = 99.;
  TLorentzVector * photon_p4 = (TLorentzVector*)pho_p4->At(photonind);
  TVector3 * photon_calopos = (TVector3*)pho_calopos->At(photonind);
  for(int iel=0;iel!=el_std_n;++iel) {
    if(el_std_hp_expin[iel]>maxlosthits)continue;
    if(el_std_scind[iel] == pho_scind[photonind]) {
      elind = iel;
      break;
    }
  }
  if(elind >= 0)
    eldr = sqrt(el_std_detain[elind]*el_std_detain[elind] + el_std_dphiin[elind]*el_std_dphiin[elind]);
  return eldr;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::IsoEcalHitsSumEtNumCrystal( TVector3 *calopos, Float_t innerConeDR, Float_t outerConeDR, Float_t stripEtaHalfWidth, Float_t stripHalfLength) {
  if(LDEBUG)std::cout << "IsoEcalHitsSumEtNumCrystal begin" << std::endl;
  Float_t ecalhitsSumEt=0.;
  for(int i=0; i!= ecalhit_n; ++i) {
    //if(ecalhit_type[i]==2)continue;// don't use preshower (still want this?? - crashes without it)
    TLorentzVector * rechitp4 = (TLorentzVector *) ecalhit_p4->At(i);
    Float_t etaclus = calopos->Eta();
    Float_t phiclus = calopos->Phi();
    Float_t eta = rechitp4->Eta();
    Float_t phi = rechitp4->Phi();
    Float_t etaDiff = eta - etaclus;
    Float_t phiDiff = phi - phiclus;
    if(phiDiff>TMath::Pi())phiDiff = TMath::TwoPi() - phiDiff;
    Float_t deltaR = calopos->DeltaR(rechitp4->Vect());
    Float_t deltaEta = fabs(rechitp4->Eta() - calopos->Eta());
    if(sqrt(etaDiff*etaDiff + phiDiff*phiDiff)>outerConeDR)continue;
    if(fabs(etaclus) < 1.479) {
      if(fabs(etaDiff) < 0.0174*stripEtaHalfWidth)continue;
      if(sqrt(etaDiff*etaDiff + phiDiff*phiDiff) < 0.0174*innerConeDR)continue;
    } else {
      if(fabs(etaDiff) < 0.00864*fabs(sinh(eta))*stripEtaHalfWidth)continue;
      if(sqrt(etaDiff*etaDiff + phiDiff*phiDiff) < 0.00864*fabs(sinh(eta))*innerConeDR)continue;
    }
    if(fabs(rechitp4->E())>0.08 && fabs(rechitp4->Et()) > 0.)ecalhitsSumEt+=rechitp4->Et();
  }
  if(LDEBUG)std::cout << "IsoEcalHitsSumEtNumCrystal end" << std::endl;
  return ecalhitsSumEt;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::WorstSumTrackPtInCone(TLorentzVector *photon_p4, Int_t returnvtxind, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {

  Int_t worstvtxind = -1;
  Float_t maxisosum = -100;
  for(int ivtx=0;ivtx!=vtx_std_n;++ivtx) {
    Float_t thisvtxisosum = SumTrackPtInCone(photon_p4, ivtx, PtMin, OuterConeRadius, InnerConeRadius, EtaStripHalfWidth, dzmax, dxymax);
    if(thisvtxisosum > maxisosum) {
      maxisosum = thisvtxisosum;
      worstvtxind = ivtx;
    }
  }

  if(returnvtxind == 1) {
    return 0.5+(float)worstvtxind;
  } else {
    return maxisosum;
  }

}


// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::SumTrackPtInCone(TLorentzVector *photon_p4, Int_t vtxind, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {
  // TRACKER Isolation
  if(vtxind<0)return -99;
  TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(vtxind);
  float SumTrackPt=0;
  for(unsigned int itk=0; itk!=tk_n; itk++) {
    TLorentzVector * tkp4= (TLorentzVector *) tk_p4->At(itk);
    if(tkp4->Pt() < PtMin)continue;
    TVector3 * tkpos= (TVector3 *) tk_vtx_pos->At(itk);
    double deltaz = fabs(vtxpos->Z() - tkpos->Z());
    if(deltaz > dzmax)continue;
    double dxy = ( -(tkpos->X() - vtxpos->X())*tkp4->Py() + (tkpos->Y() - vtxpos->Y())*tkp4->Px()) / tkp4->Pt();
    if(fabs(dxy) > dxymax)continue;
    double tk_eta = tkp4->Eta();
    double tk_phi = tkp4->Phi();
    double deta = fabs(photon_p4->Eta() - tk_eta);
    double dphi = fabs(photon_p4->Phi() - tk_phi);
    if(dphi > TMath::Pi())dphi = TMath::TwoPi() - dphi;
    double deltaR = sqrt(deta*deta + dphi*dphi);
    if(deltaR < OuterConeRadius && deltaR >= InnerConeRadius && deta >= EtaStripHalfWidth)SumTrackPt+=tkp4->Pt();
  }
  return SumTrackPt;
}


// CiC SELECTION CODE END - SSIMON

//
// Generate dictionary entries for branches from GeneralFunctions_h 
//
// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::DefineUserBranches() 
{
#ifndef __CINT__

	BRANCH_DICT(gv_n  );
	BRANCH_DICT(gv_pos);
	BRANCH_DICT(pu_n);
	BRANCH_DICT(pu_zpos);
	BRANCH_DICT(pu_sumpt_lowpt);
	BRANCH_DICT(pu_sumpt_highpt);
	BRANCH_DICT(pu_ntrks_lowpt);
	BRANCH_DICT(pu_ntrks_highpt);
	
	BRANCH_DICT(rho);
	
	BRANCH_DICT(vtx_std_sel);
	BRANCH_DICT(vtx_std_ranked_list);

	BRANCH_DICT(pho_tkiso_recvtx_030_002_0000_10_01);
	BRANCH_DICT(pho_tkiso_badvtx_040_002_0000_10_01);
	BRANCH_DICT(pho_drtotk_25_99);


	BRANCH_DICT(pho_cutlevel_lead);
	BRANCH_DICT(pho_passcuts_lead);
	BRANCH_DICT(pho_cutlevel_sublead);
	BRANCH_DICT(pho_passcuts_sublead);

#endif
}
