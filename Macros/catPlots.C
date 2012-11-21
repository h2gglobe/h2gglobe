#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TBox.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"

using namespace std;

TCanvas *canv = new TCanvas();
  
TLine l1;
TBox b1;

TH1F* makeHists(TTree *tree, string name){
  
  TH1F *total = new TH1F(Form("total_%s",name.c_str()),    name.c_str(),100,-1,1);
  TH1F *ebebhr9 = new TH1F(Form("ebebhr9_%s",name.c_str()),name.c_str(),100,-1,1);
  TH1F *ebeblr9 = new TH1F(Form("ebeblr9_%s",name.c_str()),name.c_str(),100,-1,1);
  TH1F *ebeehr9 = new TH1F(Form("ebeehr9_%s",name.c_str()),name.c_str(),100,-1,1);
  TH1F *ebeelr9 = new TH1F(Form("ebeelr9_%s",name.c_str()),name.c_str(),100,-1,1);

  float diphoton_bdt_;
  float evweight_;
  double lead_calo_eta_;
  double sublead_calo_eta_;
  float lead_r9_;
  float sublead_r9_;

  tree->SetBranchAddress("diphoton_bdt",&diphoton_bdt_);
  tree->SetBranchAddress("evweight",&evweight_);
  tree->SetBranchAddress("lead_calo_eta",&lead_calo_eta_);
  tree->SetBranchAddress("sublead_calo_eta",&sublead_calo_eta_);
  tree->SetBranchAddress("lead_r9",&lead_r9_);
  tree->SetBranchAddress("sublead_r9",&sublead_r9_);

  for (int entry=0; entry<tree->GetEntries(); entry++){
    if (entry%5000==0) cout << entry << endl;
    tree->GetEntry(entry);
    total->Fill(diphoton_bdt_,evweight_);
    // both barrel
    if (TMath::Abs(lead_calo_eta_)<1.444 && TMath::Abs(sublead_calo_eta_)<1.444){
      if (lead_r9_>0.94 && sublead_r9_>0.94) ebebhr9->Fill(diphoton_bdt_,evweight_);
      else ebeblr9->Fill(diphoton_bdt_,evweight_);
    }
    // one endcap
    else{
      if (lead_r9_>0.94 && sublead_r9_>0.94) ebeehr9->Fill(diphoton_bdt_,evweight_);
      else ebeelr9->Fill(diphoton_bdt_,evweight_);
    }
  }

  total->SetLineColor(kBlack);
  ebebhr9->SetLineColor(kRed);
  ebeblr9->SetLineColor(kGreen);
  ebeehr9->SetLineColor(kBlue);
  ebeelr9->SetLineColor(kMagenta);
  total->SetLineWidth(2);
  ebebhr9->SetLineWidth(2);
  ebeblr9->SetLineWidth(2);
  ebeehr9->SetLineWidth(2);
  ebeelr9->SetLineWidth(2);

  double ylow=0.;
  double yhigh=total->GetMaximum()*1.1;
  double x1leg=0.2;
  double x2leg=0.5;
  if (name=="data") {
    ylow=5.;
    yhigh=4e4;
    x1leg=0.5;
    x2leg=0.8;
    canv->SetLogy();
  }
  else canv->SetLogy(false);
  

  TLegend *leg = new TLegend(x1leg,0.7,x2leg,0.88);
  leg->SetFillColor(0);
  leg->AddEntry(total,"All","L");
  leg->AddEntry(ebebhr9,"both EB, both R9>0.94","L");
  leg->AddEntry(ebeblr9,"both EB, !both R9>0.94","L");
  leg->AddEntry(ebeehr9,"!both EB, both R9>0.94","L");
  leg->AddEntry(ebeelr9,"!both EB, !both R9>0.94","L");

  total->GetYaxis()->SetRangeUser(ylow,yhigh);
  total->Draw();
  ebebhr9->Draw("same");
  ebeblr9->Draw("same");
  ebeehr9->Draw("same");
  ebeelr9->Draw("same");
  
  l1.DrawLine(-0.05,ylow,-0.05,yhigh);
  l1.DrawLine(0.49,ylow,0.49,yhigh);
  l1.DrawLine(0.79,ylow,0.79,yhigh);
  l1.DrawLine(0.91,ylow,0.91,yhigh);
  b1.DrawBox(-1.,ylow,-0.05,yhigh);
  
  leg->Draw("same");
  
  canv->Print(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.C",name.c_str()));
  return total;

}


void catPlots(string infilename, bool batch=true){
  
  gROOT->SetBatch(batch);
  gStyle->SetOptStat(0);

  l1.SetLineColor(kBlue);
  l1.SetLineStyle(7);
  l1.SetLineWidth(2);
  b1.SetFillColor(kBlue);
  b1.SetFillStyle(3003);
  
  
  TFile *inFile = TFile::Open(infilename.c_str());
  
  TTree *data = (TTree*)inFile->Get("spin_trees/Data");
  TTree *ggh = (TTree*)inFile->Get("spin_trees/ggh_m125_8TeV");
  TTree *vbf = (TTree*)inFile->Get("spin_trees/vbf_m125_8TeV");
  TTree *wzh = (TTree*)inFile->Get("spin_trees/wzh_m125_8TeV");
  TTree *tth = (TTree*)inFile->Get("spin_trees/tth_m125_8TeV");
  TTree *grav = (TTree*)inFile->Get("spin_trees/rsGravGG_m125_8TeV");

  TH1F *dataH = makeHists(data,"data");
  TH1F *gghH = makeHists(ggh,"ggh");
  TH1F *vbfH = makeHists(vbf,"vbf");
  TH1F *wzhH = makeHists(wzh,"wzh");
  TH1F *tthH = makeHists(tth,"tth");
  TH1F *gravH = makeHists(grav,"grav");

  gghH->SetLineColor(kBlue);
  vbfH->SetLineColor(kRed);
  wzhH->SetLineColor(kGreen);
  tthH->SetLineColor(kMagenta);
  gravH->SetLineColor(kBlack);

  double ylow=0.;
  double yhigh=gghH->GetMaximum()*1.2;
  gghH->GetYaxis()->SetRangeUser(ylow,yhigh);
  gghH->Draw();
  vbfH->Draw("same");
  wzhH->Draw("same");
  tthH->Draw("same");
  gravH->Draw("same");
  cout << "ggh: " << gghH->GetSumOfWeights() << endl;
  cout << "vbf: " << vbfH->GetSumOfWeights() << endl;
  cout << "wzh: " << wzhH->GetSumOfWeights() << endl;
  cout << "tth: " << tthH->GetSumOfWeights() << endl;
  cout << "grav: " << gravH->GetSumOfWeights() << endl;
  
  TLegend *leg = new TLegend(0.2,0.7,0.5,0.88);
  leg->SetFillColor(0);
  leg->AddEntry(gghH,"ggh","L");
  leg->AddEntry(vbfH,"vbf","L");
  leg->AddEntry(wzhH,"wzh","L");
  leg->AddEntry(tthH,"tth","L");
  leg->AddEntry(gravH,"grav","L");
  
  l1.DrawLine(-0.05,ylow,-0.05,yhigh);
  l1.DrawLine(0.49,ylow,0.49,yhigh);
  l1.DrawLine(0.79,ylow,0.79,yhigh);
  l1.DrawLine(0.91,ylow,0.91,yhigh);
  b1.DrawBox(-1.,ylow,-0.05,yhigh);
  
  leg->Draw("same");
  
  canv->Print("all.pdf");
  canv->Print("all.C");
}
