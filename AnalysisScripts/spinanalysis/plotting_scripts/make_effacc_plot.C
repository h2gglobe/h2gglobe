#include <iostream>

float getNEvents(TTree *tree, TH1F *h){
  h->Reset();
  tree->Draw("abs(costheta_cs)>>h","evweight*(evweight>=0.)","goff");
  return h->GetSumOfWeights();
}

void make_effacc_plot(string filename, float lumi=19620., double mass=125.){
  
  gSystem->Load("../../libLoopAll.so");

  TFile *inFile = TFile::Open(filename.c_str());

  TTree *gghTree = (TTree*)inFile->Get("spin_trees/ggh_m125_8TeV");
  TTree *vbfTree = (TTree*)inFile->Get("spin_trees/vbf_m125_8TeV");
  TTree *wzhTree = (TTree*)inFile->Get("spin_trees/wzh_m125_8TeV");
  TTree *tthTree = (TTree*)inFile->Get("spin_trees/tth_m125_8TeV");

  TTree *ggh_gravTree = (TTree*)inFile->Get("spin_trees/ggh_grav_m125_8TeV");
  TTree *vbf_gravTree = (TTree*)inFile->Get("spin_trees/vbf_grav_m125_8TeV");
  
  TH1F *h = new TH1F("h","",1,0,1);
  
  float gghEvs = getNEvents(gghTree,h);
  float vbfEvs = getNEvents(vbfTree,h);
  float wzhEvs = getNEvents(wzhTree,h);
  float tthEvs = getNEvents(tthTree,h);
  float ggh_gravEvs = getNEvents(ggh_gravTree,h);
  float vbf_gravEvs = getNEvents(vbf_gravTree,h);

  float smEvs = gghEvs+vbfEvs+wzhEvs+tthEvs;

  Normalization_8TeV *norm = new Normalization_8TeV();
  float gghNorm = norm->GetBR(mass)*norm->GetXsection(mass,"ggh")*lumi;
  float vbfNorm = norm->GetBR(mass)*norm->GetXsection(mass,"vbf")*lumi;
  float wzhNorm = norm->GetBR(mass)*norm->GetXsection(mass,"wzh")*lumi;
  float tthNorm = norm->GetBR(mass)*norm->GetXsection(mass,"tth")*lumi;
  float ggh_gravNorm = norm->GetBR(mass)*norm->GetXsection(mass,"ggh_grav")*lumi;
  float vbf_gravNorm = norm->GetBR(mass)*norm->GetXsection(mass,"vbf_grav")*lumi;
  
  cout << "ggh: " << gghEvs << " " << gghNorm << endl;
  cout << "vbf: " << vbfEvs << " " << vbfNorm << endl;
  cout << "wzh: " << wzhEvs << " " << wzhNorm << endl;
  cout << "tth: " << tthEvs << " " << tthNorm << endl;
  cout << "ggh_grav: " << ggh_gravEvs << " " << ggh_gravNorm << endl;
  cout << "vbf_grav: " << vbf_gravEvs << " " << vbf_gravNorm << endl;

  TH1F *gghEffAcc = new TH1F("gghEffAcc","gghEffAcc",20,0.,1.);
  TH1F *ggh_gravEffAcc = new TH1F("ggh_gravEffAcc","ggh_gravEffAcc",20,0.,1.);

  for (int b=1; b<=gghEffAcc->GetNbinsX(); b++){
    float low = gghEffAcc->GetBinLowEdge(b);
    float high = gghEffAcc->GetBinLowEdge(b+1);
    h->Reset();
    gghTree->Draw("abs(costheta_cs)>>h",Form("evweight*((costheta_cs>=%1.3f&&costheta_cs<%1.3f))",low,high),"goff");
    gghEffAcc->SetBinContent(b,h->GetSumOfWeights()/gghNorm);
    h->Reset();
    ggh_gravTree->Draw("abs(costheta_cs)>>h",Form("evweight*((costheta_cs>=%1.3f&&costheta_cs<%1.3f))",low,high),"goff");
    ggh_gravEffAcc->SetBinContent(b,h->GetSumOfWeights()/ggh_gravNorm);
  }

  gghEffAcc->SetLineColor(kRed);
  gghEffAcc->SetLineWidth(2);
  ggh_gravEffAcc->SetLineColor(kBlue);
  ggh_gravEffAcc->SetLineWidth(2);

  gghEffAcc->Draw("HIST");
  ggh_gravEffAcc->Draw("HISTsame");


}
