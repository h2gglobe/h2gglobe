#include "TFile.h"
#include "TTree.h"
#include "TString.h"

int PhotonIDCategory(float r9, float eta) {
  //std::cout << "CAT" << std::endl;
  int thisCat = -1;
  
  bool etaCat = fabs(eta)>1.479;
  Int_t r9Cat  = (Int_t)(r9<0.90) + (Int_t)(r9<0.94);

  thisCat = 3*(etaCat) + r9Cat;
  
  return thisCat;
}

converter_ID(bool leading = false, bool subleading = false) {
  gROOT->Reset(); 

  Int_t itype;
  Float_t rho, w, nvtx;
  Float_t leadr9, leadsieie, leadhcaldr04, leadecaldr03, leadecaldr04, leadeta, leadpt;
  Float_t lead_tkiso_recvtx_030_002_0000_10_01, lead_tkiso_badvtx_040_002_0000_10_01, lead_drtotk_25_99, leadhovere;
  Float_t subleadr9, subleadsieie, subleadhcaldr04, subleadecaldr03, subleadecaldr04, subleadeta, subleadpt;
  Float_t sublead_tkiso_recvtx_030_002_0000_10_01, sublead_tkiso_badvtx_040_002_0000_10_01, sublead_drtotk_25_99, subleadhovere;
  Float_t subleadpfisocharged03, subleadpfisophoton03, subleadpfisoneutral03;
  Float_t subleadpfisochargedbadvtx04, subleadpfisophotonbadvtx04, subleadpfisoneutralbadvtx04;
  Float_t leadpfisocharged03, leadpfisophoton03, leadpfisoneutral03;
  Float_t leadpfisochargedbadvtx04, leadpfisophotonbadvtx04, leadpfisoneutralbadvtx04;
  Float_t leadcic6pfindex, subleadcic6pfindex, mass;
  Float_t scwidtheta_l, scwidthphi_l, scwidtheta_sl, scwidthphi_sl;

  Int_t cat;
  Float_t sieie, r9, goodpf_iso, badpf_iso, cic6pfindex, eta, isLeading;
  Float_t good_iso, bad_iso, drtotk, hoe, tkiso, tkisopf, pt, ptom;
  Float_t scwidtheta, scwidthphi;

  TString SEToMLOW = "0.25";
  TString LEToMLOW = "0.333";

  TString presel = "leadpt/mass>"+LEToMLOW+"&&subleadpt/mass>"+SEToMLOW+"&&abs(subleadeta)<2.5&&abs(leadeta)<2.5";
  //  presel  = presel + "&&abs(mass-120)<35&&((itype>0&&leadgenmatch&&subleadgenmatch&&abs(mass-120)<10)||((itype<0&&leadci6cindex>3)))";
  //presel  = presel + "&&abs(mass-120)<35&&((itype>0&&leadgenmatch&&subleadgenmatch&&abs(mass-120)<10)||(itype<0))";
  presel  = presel + "&&abs(mass-120)<35&&((itype>0&&leadgenmatch&&subleadgenmatch&&abs(deltaM)<10)||(itype<0))";// && !leadgenmatch&& !subleadgenmatch))";

  //TFile *f = new TFile("sept15_VBF_preplots_mutightisoPF.root");
  //TFile *f = new TFile("nov3_even_more_signal_deltaM.root");
  TFile *f = new TFile("jan14.root");
  //dec1.root");
  TTree *tinput = (TTree*)f->Get("ntuple");
  
  std::cout << "Selecting events..." << std::endl;

  TFile* tempoutput = new TFile("temp_pres.root", "RECREATE");
  TTree *temp = tinput->CopyTree(presel);
  temp->Write();
  temp->SetBranchAddress("itype", &itype);
  temp->SetBranchAddress("rho", &rho);
  temp->SetBranchAddress("nvtx", &nvtx);
  temp->SetBranchAddress("w", &w);
  temp->SetBranchAddress("mass", &mass);
  temp->SetBranchAddress("leadr9", &leadr9);
  temp->SetBranchAddress("leadsieie", &leadsieie);
  temp->SetBranchAddress("leadhcaldr04", &leadhcaldr04);
  temp->SetBranchAddress("leadecaldr03", &leadecaldr03);
  temp->SetBranchAddress("leadecaldr04", &leadecaldr04);
  temp->SetBranchAddress("leadeta", &leadeta);
  temp->SetBranchAddress("leadpt", &leadpt);
  temp->SetBranchAddress("leadci6cpfindex", &leadcic6pfindex);
  //temp->SetBranchAddress("lead_tkiso_recvtx_030_002_0000_10_01", &lead_tkiso_recvtx_030_002_0000_10_01);
  //temp->SetBranchAddress("lead_tkiso_badvtx_040_002_0000_10_01", &lead_tkiso_badvtx_040_002_0000_10_01);
  temp->SetBranchAddress("lead_drtotk_25_99", &lead_drtotk_25_99);
  temp->SetBranchAddress("lead_pfiso_charged03", &leadpfisocharged03);
  temp->SetBranchAddress("lead_pfiso_photon03", &leadpfisophoton03);
  temp->SetBranchAddress("lead_pfiso_neutral03", &leadpfisoneutral03);
  temp->SetBranchAddress("lead_pfiso_charged_badvtx_04", &leadpfisochargedbadvtx04);
  temp->SetBranchAddress("lead_pfiso_photon_badvtx_04", &leadpfisophotonbadvtx04);
  temp->SetBranchAddress("lead_pfiso_neutral_badvtx_04", &leadpfisoneutralbadvtx04);
  temp->SetBranchAddress("leadhovere", &leadhovere);
  temp->SetBranchAddress("subleadr9", &subleadr9);
  temp->SetBranchAddress("subleadsieie", &subleadsieie);
  temp->SetBranchAddress("subleadhcaldr04", &subleadhcaldr04);
  temp->SetBranchAddress("subleadecaldr03", &subleadecaldr03);
  temp->SetBranchAddress("subleadecaldr04", &subleadecaldr04);
  temp->SetBranchAddress("subleadeta", &subleadeta);
  temp->SetBranchAddress("subleadpt", &subleadpt);
  temp->SetBranchAddress("subleadci6cpfindex", &subleadcic6pfindex);
  //temp->SetBranchAddress("sublead_tkiso_recvtx_030_002_0000_10_01", &sublead_tkiso_recvtx_030_002_0000_10_01);
  //temp->SetBranchAddress("sublead_tkiso_badvtx_040_002_0000_10_01", &sublead_tkiso_badvtx_040_002_0000_10_01);
  temp->SetBranchAddress("sublead_pfiso_charged03", &subleadpfisocharged03);
  temp->SetBranchAddress("sublead_pfiso_photon03", &subleadpfisophoton03);
  temp->SetBranchAddress("sublead_pfiso_neutral03", &subleadpfisoneutral03);
  temp->SetBranchAddress("sublead_pfiso_charged_badvtx_04", &subleadpfisochargedbadvtx04);
  temp->SetBranchAddress("sublead_pfiso_photon_badvtx_04", &subleadpfisophotonbadvtx04);
  temp->SetBranchAddress("sublead_pfiso_neutral_badvtx_04", &subleadpfisoneutralbadvtx04);
  temp->SetBranchAddress("sublead_drtotk_25_99", &sublead_drtotk_25_99);
  temp->SetBranchAddress("subleadhovere", &subleadhovere);
  temp->SetBranchAddress("scwidtheta_l", &scwidtheta_l);
  temp->SetBranchAddress("scwidthphi_l", &scwidthphi_l);
  temp->SetBranchAddress("scwidtheta_sl", &scwidtheta_sl);
  temp->SetBranchAddress("scwidthphi_sl", &scwidthphi_sl);

  float rhofacpf[6]={0.075, 0.082, 0.143, 0.050, 0.091, 0.106};
  float rhofacbadpf[6]={0.141, 0.149, 0.208, 0.135, 0.162, 0.165};

  TFile* output = new TFile("prova.root","RECREATE");
  TTree *toutput = new TTree("events", "events");
  toutput->Branch("itype",       &itype,       "itype/I");
  toutput->Branch("isLeading",   &isLeading,   "isLeading/F");
  toutput->Branch("w",           &w,           "w/F");
  toutput->Branch("nvtx",           &nvtx,           "nvtx/F");
  toutput->Branch("cat",         &cat,         "cat/I");
  toutput->Branch("sieie",       &sieie,       "sieie/F");
  toutput->Branch("good_iso",    &good_iso,    "good_iso/F");
  toutput->Branch("bad_iso",     &bad_iso,     "bad_iso/F");
  toutput->Branch("goodpf_iso",  &goodpf_iso,  "goodpf_iso/F");
  toutput->Branch("badpf_iso",   &badpf_iso,   "badpf_iso/F");
  toutput->Branch("drtotk",      &drtotk,      "drtotk/F");
  toutput->Branch("hoe",         &hoe,         "hoe/F");
  toutput->Branch("tkiso",       &tkiso,       "tkiso/F");
  toutput->Branch("tkisopf",     &tkisopf,     "tkisopf/F");
  toutput->Branch("r9",          &r9,          "r9/F");
  toutput->Branch("pt",          &pt,          "pt/F");
  toutput->Branch("ptom",        &ptom,        "ptom/F");
  toutput->Branch("eta",         &eta,         "eta/F");
  toutput->Branch("cic6pfindex", &cic6pfindex, "cic6pfindex/F");
  toutput->Branch("scwidtheta", &scwidtheta, "scwidtheta/F");
  toutput->Branch("scwidthphi", &scwidthphi, "scwidthphi/F");

  std::cout << "Reducing variables..." << std::endl;

  int prescale = 0;
 
  for(int i=0; i<temp->GetEntries(); i++) {
    temp->GetEntry(i);
    if (itype == 0)
      continue;
    //prescale++;
    //if (prescale == 5)
    //  prescale = 0;

    //if (prescale != 0 && itype > 0)
    //  continue;

    if (subleading) {
      if (leadcic6pfindex > 3) {
	cat = PhotonIDCategory(subleadr9, subleadeta);
	good_iso = 0;//(sublead_tkiso_recvtx_030_002_0000_10_01 + subleadecaldr03 + subleadhcaldr04 - rho*0.17)*50./subleadpt;
	bad_iso = 0;//(sublead_tkiso_badvtx_040_002_0000_10_01 + subleadecaldr04 + subleadhcaldr04 - rho*0.52)*50./subleadpt;
	tkiso = 0;//sublead_tkiso_recvtx_030_002_0000_10_01;
	
	//goodpf_iso = (subleadpfisocharged03 + subleadpfisophoton03 - rho*rhofacpf[cat]+2.8)*50./subleadpt;
	goodpf_iso = (subleadpfisocharged03 + subleadpfisophoton03 - rho*rhofacpf[cat])*50./subleadpt;
	//if (goodpf_iso > 100) 
	//  goodpf_iso = 100;
	//goodpf_iso = sqrt(fabs(goodpf_iso));

	//badpf_iso = (subleadpfisochargedbadvtx04 + subleadpfisophotonbadvtx04 - rho*rhofacbadpf[cat]+4.8)*50./subleadpt;
	badpf_iso = (subleadpfisochargedbadvtx04 + subleadpfisophotonbadvtx04 - rho*rhofacbadpf[cat])*50./subleadpt;
	//if (badpf_iso > 100) 
	//  badpf_iso = 100;
	//badpf_iso = sqrt(fabs(badpf_iso));

	tkisopf = subleadpfisocharged03*50./subleadpt;
	//if (tkisopf > 100) 
	//  tkisopf = 100;
	//tkisopf = sqrt(fabs(tkisopf));

	drtotk = sublead_drtotk_25_99;
	//if (drtotk == 99)
	//  drtotk = 1;

	hoe = subleadhovere;
	pt = subleadpt;
	eta = fabs(subleadeta);
	sieie = subleadsieie;
	r9 = subleadr9;
	if (itype != 0) {
	  r9 = 1.0035*r9;
	  if (eta < 1.479) 
	    sieie = 0.87*sieie + 0.0011;
	  else
	    sieie = 0.99*sieie;
	}
	ptom = subleadpt/mass;
	cic6pfindex = subleadcic6pfindex;
	scwidtheta = scwidtheta_sl;
	scwidthphi = scwidthphi_sl;
	isLeading = -1.;
	toutput->Fill();
      }
    }

    if (leading) {
      if (subleadcic6pfindex > 3) {
	cat = PhotonIDCategory(leadr9, leadeta);
	good_iso = 0;//(lead_tkiso_recvtx_030_002_0000_10_01 + leadecaldr03 + leadhcaldr04 - rho*0.17)*50./leadpt;
	bad_iso = 0;//(lead_tkiso_badvtx_040_002_0000_10_01 + leadecaldr04 + leadhcaldr04 - rho*0.52)*50./leadpt;
	tkiso = 0;//lead_tkiso_recvtx_030_002_0000_10_01;

	//goodpf_iso = (leadpfisocharged03 + leadpfisophoton03 - rho*rhofacpf[cat]+2.8)*50./leadpt;
	goodpf_iso = (leadpfisocharged03 + leadpfisophoton03 - rho*rhofacpf[cat])*50./leadpt;
	//if (goodpf_iso > 100) 
	//  goodpf_iso = 100;
	//goodpf_iso = sqrt(fabs(goodpf_iso));

	//badpf_iso = (leadpfisochargedbadvtx04 + leadpfisophotonbadvtx04 - rho*rhofacbadpf[cat]+4.8)*50./leadpt;
	badpf_iso = (leadpfisochargedbadvtx04 + leadpfisophotonbadvtx04 - rho*rhofacbadpf[cat])*50./leadpt;
	//if (badpf_iso > 100) 
	//  badpf_iso = 100;
	//badpf_iso = sqrt(fabs(badpf_iso));
	
	tkisopf = leadpfisocharged03*50./leadpt;	
	//if (tkisopf > 100) 
	//  tkisopf = 100;
	//tkisopf = sqrt(fabs(tkisopf));
	
	drtotk = lead_drtotk_25_99;	
	//if (drtotk == 99)
	//  drtotk = 1;
	hoe = leadhovere;

	sieie = leadsieie;
	pt = leadpt;
	eta = fabs(leadeta);
	r9 = leadr9;
	if (itype != 0) {
	  r9 = 1.0035*r9;
	  if (eta < 1.479) 
	    sieie = 0.87*sieie + 0.0011;
	  else
	    sieie = 0.99*sieie;
	}
	ptom = leadpt/mass;
	cic6pfindex = leadcic6pfindex;
	scwidtheta = scwidtheta_l;
	scwidthphi = scwidthphi_l;
	isLeading = 1.;
	toutput->Fill();
      }
    }
  }
  
  toutput->Write();
  
  delete temp;
  delete tinput;

  tempoutput->Close();
  f->Close();

  std::cout << "Producing data, signal and bg trees..." << std::endl;

  // DATA
  //TFile* fdata = new TFile("data_ID.root", "RECREATE");
  //TTree* tdata = toutput->CopyTree("itype == 0");
  //tdata->Write();
  //fdata->Close();

  // SIGNAL
  TFile* fdata = new TFile("signal_ID.root", "RECREATE");
  TTree* tdata = toutput->CopyTree("itype > 0");
  tdata->Write();
  fdata->Close();

  // DATA
  TFile* fdata = new TFile("bg_ID.root", "RECREATE");
  TTree* tdata = toutput->CopyTree("itype < 0");
  tdata->Write();
  fdata->Close();
  
  output->Close();
}
