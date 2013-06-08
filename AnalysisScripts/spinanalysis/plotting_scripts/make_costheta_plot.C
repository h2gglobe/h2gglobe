#include <iostream>

void make_costheta_plot(string beforeCutsFile, string afterCutsFile){

  // hideous copy-paste coding - yuck!
  // run make_mc_truth_tree.C first

  float cosThetaBoundaries[4] = {0.2,0.375,0.55,0.75};

  TFile *outFile = new TFile("CosThetaOut.root","RECREATE");

  gStyle->SetOptStat(0);
  TFile *beforeCuts = TFile::Open(beforeCutsFile.c_str());
  
  TFile *afterCuts = TFile::Open(afterCutsFile.c_str());

  TTree *beforeTree = (TTree*)beforeCuts->Get("CosThetaTree");

  TTree *gghTree = (TTree*)afterCuts->Get("spin_trees/ggh_m125_8TeV");
  TTree *vbfTree = (TTree*)afterCuts->Get("spin_trees/vbf_m125_8TeV");
  TTree *wzhTree = (TTree*)afterCuts->Get("spin_trees/wzh_m125_8TeV");
  TTree *tthTree = (TTree*)afterCuts->Get("spin_trees/tth_m125_8TeV");
  //TTree *ggh_gravTree = (TTree*)afterCuts->Get("spin_trees/ggh_grav_m125_8TeV");
  //TTree *vbf_gravTree = (TTree*)afterCuts->Get("spin_trees/vbf_grav_m125_8TeV");
  TTree *ggh_gravTree = (TTree*)afterCuts->Get("spin_trees/gg_grav_m125_8TeV");
  TTree *vbf_gravTree = (TTree*)afterCuts->Get("spin_trees/qq_grav_m125_8TeV");

  TH1F *gghHistBefore = new TH1F("gghHistBefore","",40,0,1);
  TH1F *vbfHistBefore = new TH1F("vbfHistBefore","",40,0,1);
  TH1F *ggh_gravHistBefore = new TH1F("ggh_gravHistBefore","",40,0,1);
  TH1F *vbf_gravHistBefore = new TH1F("vbf_gravHistBefore","",40,0,1);

  TH1F *gghHistAfter = new TH1F("gghHistAfter","",40,0,1);
  TH1F *vbfHistAfter = new TH1F("vbfHistAfter","",40,0,1);
  TH1F *ggh_gravHistAfter = new TH1F("ggh_gravHistAfter","",40,0,1);
  TH1F *vbf_gravHistAfter = new TH1F("vbf_gravHistAfter","",40,0,1);

  TH1F *gghHistEffAcc = new TH1F("gghHistEffAcc","",40,0,1);
  TH1F *vbfHistEffAcc = new TH1F("vbfHistEffAcc","",40,0,1);
  TH1F *ggh_gravHistEffAcc = new TH1F("ggh_gravHistEffAcc","",40,0,1);
  TH1F *vbf_gravHistEffAcc = new TH1F("vbf_gravHistEffAcc","",40,0,1);
  
  gghHistBefore->Sumw2();
  vbfHistBefore->Sumw2();
  ggh_gravHistBefore->Sumw2();
  vbf_gravHistBefore->Sumw2();

  gghHistAfter->Sumw2();
  vbfHistAfter->Sumw2();
  ggh_gravHistAfter->Sumw2();
  vbf_gravHistAfter->Sumw2();

  gghHistEffAcc->Sumw2();
  vbfHistEffAcc->Sumw2();
  ggh_gravHistEffAcc->Sumw2();
  vbf_gravHistEffAcc->Sumw2();

  beforeTree->Draw("abs(cosThetaCS)>>gghHistBefore","type==2","goff");
  beforeTree->Draw("abs(cosThetaCS)>>vbfHistBefore","type==4","goff");
  beforeTree->Draw("abs(cosThetaCS)>>ggh_gravHistBefore","type==6","goff");
  beforeTree->Draw("abs(cosThetaCS)>>vbf_gravHistBefore","type==8","goff");

  TH1F *smHistBefore = (TH1F*)gghHistBefore->Clone("smHistBefore");
  smHistBefore->Add(vbfHistBefore);

  TH1F *smHistAfter = new TH1F("smHistAfter","smHistAfter",40,0,1);
  smHistAfter->Sumw2();
  gghTree->Draw("abs(costheta_cs)>>+smHistAfter","evweight*(evweight>0.)","goff");
  vbfTree->Draw("abs(costheta_cs)>>+smHistAfter","evweight*(evweight>0.)","goff");
  wzhTree->Draw("abs(costheta_cs)>>+smHistAfter","evweight*(evweight>0.)","goff");
  tthTree->Draw("abs(costheta_cs)>>+smHistAfter","evweight*(evweight>0.)","goff");

  /*
  gghTree->Draw("abs(costheta_cs)>>gghHistAfter","evweight*(evweight>0.)","goff");
  vbfTree->Draw("abs(costheta_cs)>>vbfHistAfter","evweight*(evweight>0.)","goff");
  ggh_gravTree->Draw("abs(costheta_cs)>>ggh_gravHistAfter","evweight*(evweight>0.)","goff");
  vbf_gravTree->Draw("abs(costheta_cs)>>vbf_gravHistAfter","evweight*(evweight>0.)","goff");
  */
  
  gghTree->Draw("abs(costheta_cs)>>gghHistAfter","","goff");
  vbfTree->Draw("abs(costheta_cs)>>vbfHistAfter","","goff");
  ggh_gravTree->Draw("abs(costheta_cs)>>ggh_gravHistAfter","","goff");
  vbf_gravTree->Draw("abs(costheta_cs)>>vbf_gravHistAfter","","goff");

  for (int b=1; b<=gghHistBefore->GetNbinsX(); b++){
    gghHistEffAcc->SetBinContent(b,gghHistAfter->GetBinContent(b)/gghHistBefore->GetBinContent(b));
    vbfHistEffAcc->SetBinContent(b,vbfHistAfter->GetBinContent(b)/vbfHistBefore->GetBinContent(b));
    ggh_gravHistEffAcc->SetBinContent(b,ggh_gravHistAfter->GetBinContent(b)/ggh_gravHistBefore->GetBinContent(b));
    vbf_gravHistEffAcc->SetBinContent(b,vbf_gravHistAfter->GetBinContent(b)/vbf_gravHistBefore->GetBinContent(b));
  }
  
  // per cat effs
  TH1F *ratioCat0;
  TH1F *ratioCat1;
  TH1F *ratioCat2;
  TH1F *ratioCat3;

  for (int cat=0; cat<4; cat++){
    
    TCut catcut;
    if (cat==0) catcut = "(category>=0 && category<5)";
    if (cat==1) catcut = "(category>=5 && category<10)";
    if (cat==2) catcut = "(category>=10 && category<15)";
    if (cat==3) catcut = "(category>=15 && category<20)";

    TH1F *gghHistAfterCat = new TH1F(Form("gghHistAfterCat%d",cat),"",40,0,1);
    TH1F *vbfHistAfterCat = new TH1F(Form("vbfHistAfterCat%d",cat),"",40,0,1);
    TH1F *ggh_gravHistAfterCat = new TH1F(Form("ggh_gravHistAfterCat%d",cat),"",40,0,1);
    TH1F *vbf_gravHistAfterCat = new TH1F(Form("vbf_gravHistAfterCat%d",cat),"",40,0,1);
    TH1F *smHistAfterCat = new TH1F(Form("smHistAfterCat%d",cat),"",40,0,1);
    
    gghHistAfterCat->Sumw2();
    vbfHistAfterCat->Sumw2();
    ggh_gravHistAfterCat->Sumw2();
    vbf_gravHistAfterCat->Sumw2();
    smHistAfterCat->Sumw2();

    gghTree->Draw(Form("abs(costheta_cs)>>gghHistAfterCat%d",cat),catcut,"goff");
    vbfTree->Draw(Form("abs(costheta_cs)>>vbfHistAfterCat%d",cat),catcut,"goff");
    ggh_gravTree->Draw(Form("abs(costheta_cs)>>ggh_gravHistAfterCat%d",cat),catcut,"goff");
    vbf_gravTree->Draw(Form("abs(costheta_cs)>>vbf_gravHistAfterCat%d",cat),catcut,"goff");
    gghTree->Draw(Form("abs(costheta_cs)>>+smHistAfterCat%d",cat),"evweight"*catcut,"goff");
    vbfTree->Draw(Form("abs(costheta_cs)>>+smHistAfterCat%d",cat),"evweight"*catcut,"goff");
    wzhTree->Draw(Form("abs(costheta_cs)>>+smHistAfterCat%d",cat),"evweight"*catcut,"goff");
    tthTree->Draw(Form("abs(costheta_cs)>>+smHistAfterCat%d",cat),"evweight"*catcut,"goff");

    gghHistAfterCat->Scale(1./gghHistAfterCat->Integral());
    smHistAfterCat->Scale(1./smHistAfterCat->Integral());
    gghHistBefore->Scale(1./gghHistBefore->Integral());
    smHistBefore->Scale(1./smHistBefore->Integral());

    gghHistAfterCat->Divide(gghHistBefore);
    vbfHistAfterCat->Divide(vbfHistBefore);
    ggh_gravHistAfterCat->Divide(ggh_gravHistBefore);
    vbf_gravHistAfterCat->Divide(vbf_gravHistBefore);

    smHistAfterCat->Divide(smHistBefore);

    gghHistAfterCat->SetLineColor(kRed);
    gghHistAfterCat->SetMarkerColor(kRed);
    gghHistAfterCat->SetMarkerStyle(kFullCircle);
    gghHistAfterCat->SetMarkerSize(0.8);
    vbfHistAfterCat->SetLineColor(kBlue+2);
    vbfHistAfterCat->SetLineWidth(2);
    ggh_gravHistAfterCat->SetLineColor(kGreen+1);
    ggh_gravHistAfterCat->SetMarkerColor(kGreen+1);
    ggh_gravHistAfterCat->SetMarkerStyle(kFullSquare);
    ggh_gravHistAfterCat->SetMarkerSize(0.8);
    vbf_gravHistAfterCat->SetLineColor(kBlue-7);
    vbf_gravHistAfterCat->SetMarkerColor(kBlue-7);
    vbf_gravHistAfterCat->SetMarkerStyle(kFullTriangleUp);
    vbf_gravHistAfterCat->SetMarkerSize(0.9);
    
    if (cat==0) {
      ratioCat0 = (TH1F*)ggh_gravHistAfterCat->Clone(Form("ratioCat%d",cat));
      ratioCat0->Divide(smHistAfterCat);
    }
    if (cat==1) {
      ratioCat1 = (TH1F*)ggh_gravHistAfterCat->Clone(Form("ratioCat%d",cat));
      ratioCat1->Divide(smHistAfterCat);
    }
    if (cat==2) {
      ratioCat2 = (TH1F*)ggh_gravHistAfterCat->Clone(Form("ratioCat%d",cat));
      ratioCat2->Divide(smHistAfterCat);
    }
    if (cat==3) {
      ratioCat3 = (TH1F*)ggh_gravHistAfterCat->Clone(Form("ratioCat%d",cat));
      ratioCat3->Divide(smHistAfterCat);
    }
  }
  
  ratioCat0->SetLineColor(kRed+1);
  ratioCat1->SetLineColor(kAzure-2);
  ratioCat2->SetLineColor(kGreen+2);
  ratioCat3->SetLineColor(kGray+2);

  ratioCat0->SetMarkerColor(kRed+1);
  ratioCat1->SetMarkerColor(kAzure-2);
  ratioCat2->SetMarkerColor(kGreen+2);
  ratioCat3->SetMarkerColor(kGray+2);

  ratioCat0->SetMarkerStyle(kFullCircle);
  ratioCat1->SetMarkerStyle(kFullCircle);
  ratioCat2->SetMarkerStyle(kFullCircle);
  ratioCat3->SetMarkerStyle(kFullCircle);

  ratioCat0->SetLineWidth(2);
  ratioCat1->SetLineWidth(2);
  ratioCat2->SetLineWidth(2);
  ratioCat3->SetLineWidth(2);

  ratioCat0->Rebin(2);
  ratioCat1->Rebin(2);
  ratioCat2->Rebin(2);
  ratioCat3->Rebin(2);

  ratioCat0->Scale(0.5);
  ratioCat1->Scale(0.5);
  ratioCat2->Scale(0.5);
  ratioCat3->Scale(0.5);

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextSize(0.04);

  TLegend *catLeg = new TLegend(0.75,0.65,0.89,0.89);
  catLeg->SetFillColor(0);
  catLeg->SetLineColor(0);
  catLeg->AddEntry(ratioCat0,"cat0","LEP");
  catLeg->AddEntry(ratioCat1,"cat1","LEP");
  catLeg->AddEntry(ratioCat2,"cat2","LEP");
  catLeg->AddEntry(ratioCat3,"cat3","LEP");

  TCanvas *cats = new TCanvas("cats","cats",10,10,700,500);
  ratioCat0->GetYaxis()->SetRangeUser(0.,2.);
  ratioCat0->GetYaxis()->SetTitle("Acc. #times Eff ratio to SM");
  ratioCat0->GetXaxis()->SetTitle("|cos(#theta_{CS})|");
  TF1 *f = new TF1("f1","1.",0.,1.);
  f->SetLineStyle(kDashed);
  f->SetLineWidth(2);
  f->SetLineColor(kBlack);
  ratioCat0->Draw("AXIS");
  f->Draw("Lsame");
  for (int b=0; b<4; b++){
    TLine *line = new TLine(cosThetaBoundaries[b],0.,cosThetaBoundaries[b],2.);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    line->Draw("same");
  }
  ratioCat0->Draw("LEPsame");
  ratioCat1->Draw("LEPsame");
  ratioCat2->Draw("LEPsame");
  ratioCat3->Draw("LEPsame");
  catLeg->Draw("same");
  lat->DrawLatex(0.11,0.92,"CMS Preliminary");
  cats->SetTicks();
  cats->Update();
  cats->Print("cats.pdf");
 

  gghHistBefore->Scale(1./gghHistBefore->Integral());
  vbfHistBefore->Scale(1./vbfHistBefore->Integral());
  ggh_gravHistBefore->Scale(1./ggh_gravHistBefore->Integral());
  vbf_gravHistBefore->Scale(1./vbf_gravHistBefore->Integral());
  smHistBefore->Scale(1./smHistBefore->Integral());

  gghHistAfter->Scale(1./gghHistAfter->Integral());
  vbfHistAfter->Scale(1./vbfHistAfter->Integral());
  ggh_gravHistAfter->Scale(1./ggh_gravHistAfter->Integral());
  vbf_gravHistAfter->Scale(1./vbf_gravHistAfter->Integral());
  smHistAfter->Scale(1./smHistAfter->Integral());

  gghHistBefore->SetLineColor(kRed);
  gghHistBefore->SetMarkerColor(kRed);
  gghHistBefore->SetMarkerStyle(kFullCircle);
  gghHistBefore->SetMarkerSize(0.8);
  vbfHistBefore->SetLineColor(kBlue+2);
  vbfHistBefore->SetLineWidth(2);
  ggh_gravHistBefore->SetLineColor(kGreen+1);
  ggh_gravHistBefore->SetMarkerColor(kGreen+1);
  ggh_gravHistBefore->SetMarkerStyle(kFullSquare);
  ggh_gravHistBefore->SetMarkerSize(0.8);
  vbf_gravHistBefore->SetLineColor(kBlue-7);
  vbf_gravHistBefore->SetMarkerColor(kBlue-7);
  vbf_gravHistBefore->SetMarkerStyle(kFullTriangleUp);
  vbf_gravHistBefore->SetMarkerSize(0.9);
  smHistBefore->SetLineColor(kRed);
  smHistBefore->SetMarkerColor(kRed);
  smHistBefore->SetMarkerStyle(kFullCircle);
  smHistBefore->SetMarkerSize(0.8);

  gghHistAfter->SetLineColor(kRed);
  gghHistAfter->SetMarkerColor(kRed);
  gghHistAfter->SetMarkerStyle(kFullCircle);
  gghHistAfter->SetMarkerSize(0.8);
  vbfHistAfter->SetLineColor(kBlue+2);
  vbfHistAfter->SetLineWidth(2);
  ggh_gravHistAfter->SetLineColor(kGreen+1);
  ggh_gravHistAfter->SetMarkerColor(kGreen+1);
  ggh_gravHistAfter->SetMarkerStyle(kFullSquare);
  ggh_gravHistAfter->SetMarkerSize(0.8);
  vbf_gravHistAfter->SetLineColor(kBlue-7);
  vbf_gravHistAfter->SetMarkerColor(kBlue-7);
  vbf_gravHistAfter->SetMarkerStyle(kFullTriangleUp);
  vbf_gravHistAfter->SetMarkerSize(0.9);
  smHistAfter->SetLineColor(kRed);
  smHistAfter->SetMarkerColor(kRed);
  smHistAfter->SetMarkerStyle(kFullCircle);
  smHistAfter->SetMarkerSize(0.8);

  gghHistEffAcc->SetLineColor(kRed);
  gghHistEffAcc->SetMarkerColor(kRed);
  gghHistEffAcc->SetMarkerStyle(kFullCircle);
  gghHistEffAcc->SetMarkerSize(0.8);
  vbfHistEffAcc->SetLineColor(kBlue+2);
  vbfHistEffAcc->SetLineWidth(2);
  ggh_gravHistEffAcc->SetLineColor(kGreen+1);
  ggh_gravHistEffAcc->SetMarkerColor(kGreen+1);
  ggh_gravHistEffAcc->SetMarkerStyle(kFullSquare);
  ggh_gravHistEffAcc->SetMarkerSize(0.8);
  vbf_gravHistEffAcc->SetLineColor(kBlue-7);
  vbf_gravHistEffAcc->SetMarkerColor(kBlue-7);
  vbf_gravHistEffAcc->SetMarkerStyle(kFullTriangleUp);
  vbf_gravHistEffAcc->SetMarkerSize(0.9);

  TLegend *leg = new TLegend(0.11,0.7,0.3,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  //leg->AddEntry(gghHistBefore,"gg#rightarrow0^{+}","LEP");
  //leg->AddEntry(vbfHistBefore,"qq#rightarrow0^{+}","L");
  //leg->AddEntry(ggh_gravHistBefore,"gg#rightarrow2_{m}^{+}","LEP");
  //leg->AddEntry(vbf_gravHistBefore,"q#bar{q}#rightarrow2_{m}^{+}","LEP");
  leg->AddEntry(smHistBefore,"0^{+} (SM)","LEP");
  leg->AddEntry(ggh_gravHistBefore,"2_{m}^{+} (gg)","LEP");
  leg->AddEntry(vbf_gravHistBefore,"2_{m}^{+} (q#bar{q})","LEP");

  // plot truth
  TCanvas *canvBefore = new TCanvas("before","before",110,10,700,500);
  ggh_gravHistBefore->GetYaxis()->SetRangeUser(0.,0.06);
  ggh_gravHistBefore->GetYaxis()->SetTitle("a.u.");
  ggh_gravHistBefore->GetXaxis()->SetTitle("|cos(#theta_{CS})|");
  ggh_gravHistBefore->GetYaxis()->SetTitleOffset(1.2);
  ggh_gravHistBefore->Draw("AXIS");
  for (int b=0; b<4; b++){
    TLine *line = new TLine(cosThetaBoundaries[b],0.,cosThetaBoundaries[b],0.06);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    line->Draw("same");
  }
  ggh_gravHistBefore->Draw("LEPsame");
  vbf_gravHistBefore->Draw("LEPsame");
  smHistBefore->Draw("LEPsame");
  //vbfHistBefore->Draw("HISTsame");
  //gghHistBefore->Draw("LEPsame");
  leg->Draw("same");
  lat->DrawLatex(0.11,0.92,"CMS Preliminary");
  canvBefore->SetTicks();
  canvBefore->Update();
  canvBefore->Print("beforeCuts.pdf");

  // plot reco+selection
  leg->SetX1NDC(0.69);
  leg->SetX2NDC(0.89);
  TCanvas *canvAfter = new TCanvas("after","after",210,10,700,500);
  ggh_gravHistAfter->GetYaxis()->SetRangeUser(0.,0.05);
  ggh_gravHistAfter->GetYaxis()->SetTitle("a.u.");
  ggh_gravHistAfter->GetXaxis()->SetTitle("|cos(#theta_{CS})|");
  ggh_gravHistAfter->GetYaxis()->SetTitleOffset(1.2);
  ggh_gravHistAfter->Draw("AXIS");
  for (int b=0; b<4; b++){
    TLine *line = new TLine(cosThetaBoundaries[b],0.,cosThetaBoundaries[b],0.05);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    line->Draw("same");
  }
  ggh_gravHistAfter->Draw("LEPsame");
  vbf_gravHistAfter->Draw("LEPsame");
  smHistAfter->Draw("LEPsame");
  //vbfHistAfter->Draw("HISTsame");
  //gghHistAfter->Draw("LEPsame");
  leg->Draw("same");
  lat->DrawLatex(0.11,0.92,"CMS Preliminary");
  canvAfter->SetTicks();
  canvAfter->Update();
  canvAfter->Print("afterCuts.pdf");
  outFile->cd();
  ggh_gravHistAfter->Write();
  vbf_gravHistAfter->Write();
  smHistAfter->Write();

  TCanvas *canvEffAcc = new TCanvas("effacc","effacc",310,10,700,500);
  ggh_gravHistEffAcc->GetYaxis()->SetRangeUser(0.,0.7);
  ggh_gravHistEffAcc->GetYaxis()->SetTitle("#epsilon#times#alpha");
  ggh_gravHistEffAcc->GetXaxis()->SetTitle("|cos(#theta_{CS})|");
  ggh_gravHistEffAcc->GetYaxis()->SetTitleOffset(1.2);
  ggh_gravHistEffAcc->Draw("AXIS");
  for (int b=0; b<4; b++){
    TLine *line = new TLine(cosThetaBoundaries[b],0.,cosThetaBoundaries[b],0.7);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    line->Draw("same");
  }
  ggh_gravHistEffAcc->Draw("LEPsame");
  vbf_gravHistEffAcc->Draw("LEPsame");
  vbfHistEffAcc->Draw("HISTsame");
  gghHistEffAcc->Draw("LEPsame");
  canvEffAcc->Print("fullEffAcc.pdf");
  

}
