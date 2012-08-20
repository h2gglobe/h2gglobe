#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMVA/Reader.h"

#include <iostream>

Int_t itype,  diphocat2r92eta;
Float_t lpt, slpt, dppt, mass, leta, sleta, lr9, slr9, w;
Float_t sublead_pfiso_charged03, sublead_pfiso_photon03, rho, sublead_pfiso_badvtx_charged04, sublead_pfiso_badvtx_photon04, sublead_drtotk_25_99;
Float_t subleadhovere, subleadsieie, lead_pfiso_charged03, lead_pfiso_photon03, lead_pfiso_badvtx_charged04, lead_pfiso_badvtx_photon04, lead_drtotk_25_99;
Float_t scwidtheta_l, scwidthphi_l, scwidtheta_sl, scwidthphi_sl;
Float_t leadhovere, leadsieie;
Float_t etamax, dmom, dmom0, subleadci6cpfmvaptom, leadci6cpfmvaptom;
Float_t pvtx, nvtx, sigma_mz, subleadptomass, dpptom, sumptom, leadptomass;
Float_t subleadmva, leadmva;
Float_t subleadmva2, leadmva2, subleadmva3, leadmva3, subleadmva4, leadmva4;
Float_t subleadmva5, leadmva5;

Float_t tmva_sieie, tmva_goodpf_iso, tmva_badpf_iso, tmva_drtotk, tmva_nvtx;
Float_t tmva_hoe, tmva_tkisopf, tmva_r9, tmva_ptom, tmva_eta, tmva_isLeading;  
Float_t tmva_scwidtheta, tmva_scwidthphi;

float rhofacpf[6]    = {0.075, 0.082, 0.143, 0.050, 0.091, 0.106};
float rhofacbadpf[6] = {0.141, 0.149, 0.208, 0.135, 0.162, 0.165};

TMVA::Reader* tmvaReader, *tmvaReader2, *tmvaReader3, *tmvaReader4;

int PhotonIDCategory(float r9, float eta) {
  int thisCat = -1;
  
  bool etaCat = fabs(eta)>1.479;
  Int_t r9Cat  = (Int_t)(r9<0.90) + (Int_t)(r9<0.94);

  thisCat = 3*(etaCat) + r9Cat;
  
  return thisCat;
}

void SetBDT() {

  tmvaReader = new TMVA::Reader("!Color:Silent");  
  tmvaReader->AddVariable("sieie", &tmva_sieie);
  tmvaReader->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  tmvaReader->AddVariable("badpf_iso", &tmva_badpf_iso);
  tmvaReader->AddVariable("drtotk", &tmva_drtotk);
  tmvaReader->AddVariable("hoe", &tmva_hoe);
  tmvaReader->AddVariable("tkisopf", &tmva_tkisopf);
  tmvaReader->AddVariable("r9", &tmva_r9);
  tmvaReader->AddVariable("ptom", &tmva_ptom);
  tmvaReader->AddVariable("eta", &tmva_eta);
  tmvaReader->BookMVA("Gradient", "weights/id_ucsd_Gradient.weights.xml");

  tmvaReader2 = new TMVA::Reader("!Color:Silent");
  tmvaReader2->AddVariable("sieie", &tmva_sieie);
  tmvaReader2->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  tmvaReader2->AddVariable("badpf_iso", &tmva_badpf_iso);
  tmvaReader2->AddVariable("drtotk", &tmva_drtotk);
  tmvaReader2->AddVariable("hoe", &tmva_hoe);
  tmvaReader2->AddVariable("tkisopf", &tmva_tkisopf);
  tmvaReader2->AddVariable("r9", &tmva_r9);
  tmvaReader2->AddVariable("ptom", &tmva_ptom);
  tmvaReader2->AddVariable("eta", &tmva_eta);
  tmvaReader2->AddVariable("nvtx", &tmva_nvtx);
  tmvaReader2->AddVariable("scwidtheta", &tmva_scwidtheta);
  tmvaReader2->AddVariable("scwidthphi", &tmva_scwidthphi);
  tmvaReader2->BookMVA("Gradient", "weights/id_ucsd_with_scwidth_Gradient.weights.xml");

  tmvaReader3 = new TMVA::Reader("!Color:Silent");
  tmvaReader3->AddVariable("sieie", &tmva_sieie);
  tmvaReader3->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  tmvaReader3->AddVariable("badpf_iso", &tmva_badpf_iso);
  tmvaReader3->AddVariable("drtotk", &tmva_drtotk);
  tmvaReader3->AddVariable("hoe", &tmva_hoe);
  tmvaReader3->AddVariable("tkisopf", &tmva_tkisopf);
  tmvaReader3->AddVariable("r9", &tmva_r9);
  tmvaReader3->AddVariable("ptom", &tmva_ptom);
  tmvaReader3->AddVariable("eta", &tmva_eta);
  tmvaReader3->AddVariable("nvtx", &tmva_nvtx);
  tmvaReader3->AddSpectator("isLeading", &tmva_isLeading);

  tmvaReader3->BookMVA("Gradient", "weights/id_ucsd_fake_Gradient.weights.xml");

  tmvaReader4 = new TMVA::Reader("!Color:Silent");
  
  tmvaReader4->AddVariable("sieie", &tmva_sieie);
  tmvaReader4->AddVariable("goodpf_iso", &tmva_goodpf_iso);
  tmvaReader4->AddVariable("badpf_iso", &tmva_badpf_iso);
  tmvaReader4->AddVariable("drtotk", &tmva_drtotk);
  tmvaReader4->AddVariable("hoe", &tmva_hoe);
  tmvaReader4->AddVariable("tkisopf", &tmva_tkisopf);
  tmvaReader4->AddVariable("r9", &tmva_r9);
  tmvaReader4->AddVariable("eta", &tmva_eta);
  tmvaReader4->AddVariable("nvtx", &tmva_nvtx);
  tmvaReader4->AddSpectator("isLeading", &tmva_isLeading);

  tmvaReader4->BookMVA("Gradient", "weights/id_ucsd_no_ptom_fake_Gradient.weights.xml");
}

Float_t BDT(int i) {
  if (i == 1)
    return tmvaReader->EvaluateMVA("Gradient");
  if (i == 2)
    return tmvaReader2->EvaluateMVA("Gradient");
  if (i == 3)
    return tmvaReader3->EvaluateMVA("Gradient");
  if (i == 4)
    return tmvaReader4->EvaluateMVA("Gradient");
  
  return -99;
}


void converter() {
  //gROOT->Reset(); 

  TString SEToMLOW = "0.25";
  TString LEToMLOW = "0.333";
  TString presel = "(leadpt/mass)>" + LEToMLOW + " && (subleadpt/mass)>" + SEToMLOW + " && abs(subleadeta)<2.5 && abs(leadeta)<2.5";
  //copysel = presel + " && abs(mass-120)<30 && ((itype>0 && leadgenmatch && subleadgenmatch) || (itype<0) || (itype==0)) && subleadcutindex>2 && leadcutindex>2";
  TString copysel = presel + " && (mass > 90) && ((itype>0 && leadgenmatch && subleadgenmatch) || (itype<0) || (itype==0)) && subleadcutindex>2 && leadcutindex>2";

  TFile *f = new TFile("jan14.root");
  TTree *tinput = (TTree*)f->Get("ntuple");
  
  SetBDT();

  std::cout << "Selecting events..." << std::endl;

  TFile* tempoutput = new TFile("temp_pres.root", "RECREATE");
  TTree *temp = tinput->CopyTree(copysel);
  temp->Write();

  temp->SetBranchStatus("*", 0);
  temp->SetBranchStatus("diphocat2r92eta",1);
  temp->SetBranchStatus("leadpt",         1);
  temp->SetBranchStatus("subleadpt",      1);
  temp->SetBranchStatus("leadeta",        1);
  temp->SetBranchStatus("subleadeta",     1);
  temp->SetBranchStatus("leadr9",         1);
  temp->SetBranchStatus("subleadr9",      1);
  temp->SetBranchStatus("diphopt",        1);
  temp->SetBranchStatus("mass",           1);
  temp->SetBranchStatus("w",              1);
  temp->SetBranchStatus("sublead_pfiso_charged03"       ,1);
  temp->SetBranchStatus("sublead_pfiso_photon03"        ,1);
  temp->SetBranchStatus("rho"                           ,1);
  temp->SetBranchStatus("sublead_pfiso_charged_badvtx_04",1);
  temp->SetBranchStatus("sublead_pfiso_photon_badvtx_04" ,1);
  temp->SetBranchStatus("sublead_drtotk_25_99"          ,1);
  temp->SetBranchStatus("subleadhovere"                 ,1);
  temp->SetBranchStatus("subleadsieie"                  ,1);
  temp->SetBranchStatus("lead_pfiso_charged03"          ,1);
  temp->SetBranchStatus("lead_pfiso_photon03"           ,1);
  temp->SetBranchStatus("lead_pfiso_charged_badvtx_04"   ,1);
  temp->SetBranchStatus("lead_pfiso_photon_badvtx_04"    ,1);
  temp->SetBranchStatus("lead_drtotk_25_99"             ,1);
  temp->SetBranchStatus("leadhovere"                    ,1);
  temp->SetBranchStatus("leadsieie"                     ,1);
  temp->SetBranchStatus("etamax",          1);
  temp->SetBranchStatus("itype",           1);
  temp->SetBranchStatus("dmom0",            1);
  temp->SetBranchStatus("dmom",             1);
  temp->SetBranchStatus("leadci6cpfmvaptom", 1);
  temp->SetBranchStatus("subleadci6cpfmvaptom", 1);
  temp->SetBranchStatus("pvtx",            1);
  temp->SetBranchStatus("nvtx",            1);
  temp->SetBranchStatus("sigma_mz",        1);
  temp->SetBranchStatus("scwidtheta_l",    1);
  temp->SetBranchStatus("scwidthphi_l",    1);
  temp->SetBranchStatus("scwidtheta_sl",   1);
  temp->SetBranchStatus("scwidthphi_sl",   1);

  temp->SetBranchAddress("diphocat2r92eta", &diphocat2r92eta);
  temp->SetBranchAddress("leadpt",          &lpt);
  temp->SetBranchAddress("subleadpt",       &slpt);
  temp->SetBranchAddress("leadeta",         &leta);
  temp->SetBranchAddress("subleadeta",      &sleta);
  temp->SetBranchAddress("leadr9",          &lr9);
  temp->SetBranchAddress("subleadr9",       &slr9);
  temp->SetBranchAddress("diphopt",         &dppt);
  temp->SetBranchAddress("mass",            &mass);
  temp->SetBranchAddress("w",               &w);
  temp->SetBranchAddress("sublead_pfiso_charged03"       ,&sublead_pfiso_charged03);
  temp->SetBranchAddress("sublead_pfiso_photon03"        ,&sublead_pfiso_photon03);
  temp->SetBranchAddress("rho"                           ,&rho);
  temp->SetBranchAddress("sublead_pfiso_charged_badvtx_04",&sublead_pfiso_badvtx_charged04);
  temp->SetBranchAddress("sublead_pfiso_photon_badvtx_04" ,&sublead_pfiso_badvtx_photon04);
  temp->SetBranchAddress("sublead_drtotk_25_99"          ,&sublead_drtotk_25_99);
  temp->SetBranchAddress("subleadhovere"                 ,&subleadhovere);
  temp->SetBranchAddress("subleadsieie"                  ,&subleadsieie);
  temp->SetBranchAddress("lead_pfiso_charged03"          ,&lead_pfiso_charged03);
  temp->SetBranchAddress("lead_pfiso_photon03"           ,&lead_pfiso_photon03);
  temp->SetBranchAddress("lead_pfiso_charged_badvtx_04"   ,&lead_pfiso_badvtx_charged04);
  temp->SetBranchAddress("lead_pfiso_photon_badvtx_04"    ,&lead_pfiso_badvtx_photon04);
  temp->SetBranchAddress("lead_drtotk_25_99"             ,&lead_drtotk_25_99);
  temp->SetBranchAddress("leadhovere"                    ,&leadhovere);
  temp->SetBranchAddress("leadsieie"                     ,&leadsieie);
  temp->SetBranchAddress("etamax",               &etamax);
  temp->SetBranchAddress("itype",                &itype);
  temp->SetBranchAddress("dmom0",                &dmom0);
  temp->SetBranchAddress("dmom",                 &dmom);
  temp->SetBranchAddress("leadci6cpfmvaptom",    &leadci6cpfmvaptom);
  temp->SetBranchAddress("subleadci6cpfmvaptom", &subleadci6cpfmvaptom);
  temp->SetBranchAddress("pvtx",                 &pvtx);
  temp->SetBranchAddress("nvtx",                 &nvtx);
  temp->SetBranchAddress("sigma_mz",             &sigma_mz);
  temp->SetBranchAddress("scwidtheta_l",         &scwidtheta_l);
  temp->SetBranchAddress("scwidthphi_l",         &scwidthphi_l);
  temp->SetBranchAddress("scwidtheta_sl",        &scwidtheta_sl);
  temp->SetBranchAddress("scwidthphi_sl",        &scwidthphi_sl);

  TFile* output = new TFile("prova.root","RECREATE");
  TTree *toutput = new TTree("events", "events");

  toutput->Branch("w",               &w,               "w/F");
  toutput->Branch("subleadptomass",  &subleadptomass,  "subleadptomass/F");
  toutput->Branch("diphoptom",       &dpptom,          "diphoptom/F");
  toutput->Branch("sumptom",         &sumptom,         "sumptom/F");
  toutput->Branch("subleadmva",      &subleadmva,      "subleadmva/F");
  toutput->Branch("leadmva",         &leadmva,         "leadmva/F");
  toutput->Branch("subleadmva2",     &subleadmva2,     "subleadmva2/F");
  toutput->Branch("leadmva2",        &leadmva2,        "leadmva2/F");
  toutput->Branch("subleadmva3",     &subleadmva,      "subleadmva3/F");
  toutput->Branch("leadmva3",        &leadmva3,        "leadmva3/F");
  toutput->Branch("subleadmva4",     &subleadmva4,     "subleadmva4/F");
  toutput->Branch("leadmva4",        &leadmva4,        "leadmva4/F");
  toutput->Branch("subleadmva5",     &subleadmva5,     "subleadmva5/F");
  toutput->Branch("leadmva5",        &leadmva5,        "leadmva5/F");
  toutput->Branch("etamax",          &etamax,          "etamax/F");
  toutput->Branch("leadeta",         &leta,            "leadeta/F");
  toutput->Branch("subleadeta",      &sleta,           "subleadeta/F");
  toutput->Branch("leadr9",          &lr9,             "leadr9/F");
  toutput->Branch("subleadr9",       &slr9,            "subleadr9/F");
  toutput->Branch("itype",           &itype,           "itype/I");
  toutput->Branch("dmom0",           &dmom0,           "dmom0/F");
  toutput->Branch("dmom",            &dmom,            "dmom/F");
  //toutput->Branch("diphocat2r92eta", &diphocat2r92eta, "diphocat2r92eta/F");
  toutput->Branch("pvtx",            &pvtx,            "pvtx/F");
  toutput->Branch("nvtx",            &nvtx,            "nvtx/F");
  toutput->Branch("sigma_mz",        &sigma_mz,        "sigma_mz/F");
  
  std::cout << "Reducing variables..." << std::endl;
  for(int i=0; i<temp->GetEntries(); i++) {
    temp->GetEntry(i);

    if (itype == 0)
      continue;
    
    if (pvtx > 1.)
      pvtx = 1.;

    if (pvtx < 0.)
      pvtx = 0.;

    //diphocat = float(diphocat2r92eta);
    Int_t cat = PhotonIDCategory(slr9, sleta);
    tmva_isLeading = -1;
    tmva_goodpf_iso = (sublead_pfiso_charged03 + sublead_pfiso_photon03 - rho*rhofacpf[cat])*50./slpt;
    tmva_badpf_iso = (sublead_pfiso_badvtx_charged04 + sublead_pfiso_badvtx_photon04 - rho*rhofacbadpf[cat])*50./slpt;
    tmva_tkisopf = sublead_pfiso_charged03*50./slpt;
    tmva_drtotk = sublead_drtotk_25_99;
    tmva_hoe = subleadhovere;
    tmva_sieie = subleadsieie;
    tmva_r9 = slr9;
    tmva_eta = fabs(sleta);
      
    if (itype != 0) {
      tmva_r9 = 1.0035*tmva_r9;
      if (tmva_eta < 1.479) 
	tmva_sieie = 0.87*tmva_sieie + 0.0011;
      else
	tmva_sieie = 0.99*tmva_sieie;
    }

    tmva_ptom = slpt/mass;
    //    std::cout << tmva_goodpf_iso << " " << tmva_badpf_iso << " " << tmva_drtotk << " " << tmva_tkisopf << " " << tmva_hoe << " " 
    //	      <<  tmva_sieie << " " << tmva_r9 << " " << tmva_eta << " " << tmva_ptom << " " << BDT(1) << std::endl;
    tmva_scwidtheta = scwidtheta_sl;
    tmva_scwidthphi = scwidthphi_sl;
    subleadmva = BDT(1);
    subleadmva2 = BDT(2);
    subleadmva3 = BDT(3);
    subleadmva4 = BDT(4);
    subleadmva5 = subleadci6cpfmvaptom;

    cat = PhotonIDCategory(lr9, leta);
    tmva_goodpf_iso = (lead_pfiso_charged03 + lead_pfiso_photon03 - rho*rhofacpf[cat])*50./lpt;
    tmva_badpf_iso = (lead_pfiso_badvtx_charged04 + lead_pfiso_badvtx_photon04 - rho*rhofacbadpf[cat])*50./lpt;
    tmva_tkisopf = lead_pfiso_charged03*50./lpt;
    tmva_drtotk = lead_drtotk_25_99;
    tmva_hoe = leadhovere;
    tmva_sieie = leadsieie;
    tmva_eta = fabs(leta);
    tmva_r9 = lr9;
    if (itype != 0) {
      tmva_r9 = 1.0035*tmva_r9;
        if (tmva_eta < 1.479) 
	tmva_sieie = 0.87*tmva_sieie + 0.0011;
      else
	tmva_sieie = 0.99*tmva_sieie;
    }

    tmva_ptom = lpt/mass;
    //std::cout << tmva_goodpf_iso << " " << tmva_badpf_iso << " " << tmva_drtotk << " " << tmva_tkisopf << " " << tmva_hoe << " " 
    //	      <<  tmva_sieie << " " << tmva_r9 << " " << tmva_eta << " " << tmva_ptom << " " << BDT(1) << std::endl;
    //std::cout << "---------" << std::endl;
    tmva_scwidtheta = scwidtheta_l;
    tmva_scwidthphi = scwidthphi_l;
    leadmva = BDT(1);
    leadmva2 = BDT(2);
    leadmva3 = BDT(3);
    leadmva4 = BDT(4);
    leadmva5 = leadci6cpfmvaptom;

    dpptom = dppt/mass;
    subleadptomass = slpt/mass; 
    leadptomass = lpt/mass;
    sumptom = leadptomass + subleadptomass;
    leta = fabs(leta);
    sleta = fabs(sleta);
    sigma_mz = fabs(sigma_mz)/mass;

    //std::cout << dmom0 << itype << std::endl;
    if (itype > 0)
      w = w*(0.01/dmom0);
    
    toutput->Fill();
  }
  
  toutput->Write();
  
  tempoutput->Close();
  f->Close();

  std::cout << "Producing data, signal and bg trees..." << std::endl;

  // DATA
  //TFile* fdata = new TFile("data.root", "RECREATE");
  //TTree* tdata = toutput->CopyTree("itype == 0");
  //tdata->Write();
  //fdata->Close();

  // SIGNAL
  fdata = new TFile("signal.root", "RECREATE");
  tdata = toutput->CopyTree("itype > 0");
  tdata->Write();
  fdata->Close();

  // BG
  fdata = new TFile("bg.root", "RECREATE");
  tdata = toutput->CopyTree("itype < 0");
  tdata->Write();
  fdata->Close();
  
  output->Close();
}
