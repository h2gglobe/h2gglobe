#include <TROOT.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TPostScript.h"
#include <TPaveStats.h>
#include "TLegend.h"
#include <iostream>
#include <TProfile.h>
#include "TGraph.h"

void EffPlotter() {
  Double_t x1 =.72, y1= 0.72,x2 = x1+0.25, y2= y1+0.20;
  TFile *histfile  = new TFile("/afs/cern.ch/user/s/swmukher/work/FutureHiggsStudy/New/CMSSW_6_1_2/src/h2gglobe/AnalysisScripts/histograms_CMS-HGG.root");

  //////// Output File /////////
  TFile* outputFile = new TFile("Comp_8_14_TeV.root","RECREATE");

  /// Eff vs Eta ///
  TH1F* num_8TeV = (TH1F*)histfile->Get("num_hist_absEta_cat0_ggh_8TeV");
  TH1F* den_8TeV = (TH1F*)histfile->Get("deno_hist_absEta_cat0_ggh_8TeV");
  TH1F reco_eff_8TeV = ((*num_8TeV)/(*den_8TeV));
  reco_eff_8TeV.SetStats(0);
  reco_eff_8TeV.SetLineColor(kBlue);
  reco_eff_8TeV.SetLineWidth(2);
  reco_eff_8TeV.GetYaxis()->SetTitle("Reconstruction Efficiency");
  reco_eff_8TeV.SetTitle("Reco Efficiency");

  TH1F* num_14TeV = (TH1F*)histfile->Get("num_hist_absEta_cat0_ggh_14TeV_PU50");
  TH1F* den_14TeV = (TH1F*)histfile->Get("deno_hist_absEta_cat0_ggh_14TeV_PU50");
  TH1F reco_eff_14TeV = ((*num_14TeV)/(*den_14TeV));
  reco_eff_14TeV.SetStats(0);
  reco_eff_14TeV.SetLineColor(kRed);
  reco_eff_14TeV.SetLineWidth(2);

  TCanvas* PhoRecoEff = new TCanvas("Pho_reco_eff","RecoEff_compare", 800, 600);
  PhoRecoEff->cd();

  reco_eff_8TeV.DrawCopy();
  reco_eff_14TeV.DrawCopy("same");
  TLegend* L1 = new TLegend(x1,y1,x2,y2);
  L1->AddEntry(&reco_eff_8TeV,  "8TeV",    "L");
  L1->AddEntry(&reco_eff_14TeV, "14TeV",   "L");
  L1->Draw();
  PhoRecoEff->Update();
  PhoRecoEff->Write();


  /// Eff vs Pt ///
  TH1F* num_pt_8TeV = (TH1F*)histfile->Get("num_hist_Pt_cat0_ggh_8TeV");
  TH1F* den_pt_8TeV = (TH1F*)histfile->Get("deno_hist_Pt_cat0_ggh_8TeV");
  TH1F reco_eff_pt_8TeV = ((*num_pt_8TeV)/(*den_pt_8TeV));
  reco_eff_pt_8TeV.SetStats(0);
  reco_eff_pt_8TeV.SetLineColor(kBlue);
  reco_eff_pt_8TeV.SetLineWidth(2);
  reco_eff_pt_8TeV.GetYaxis()->SetTitle("Reconstruction Efficiency");
  reco_eff_pt_8TeV.SetTitle("Reco Efficiency");

  TH1F* num_pt_14TeV = (TH1F*)histfile->Get("num_hist_Pt_cat0_ggh_14TeV_PU50");
  TH1F* den_pt_14TeV = (TH1F*)histfile->Get("deno_hist_Pt_cat0_ggh_14TeV_PU50");
  TH1F reco_eff_pt_14TeV = ((*num_pt_14TeV)/(*den_pt_14TeV));
  reco_eff_pt_14TeV.SetStats(0);
  reco_eff_pt_14TeV.SetLineColor(kRed);
  reco_eff_pt_14TeV.SetLineWidth(2);

  TCanvas* PhoRecoEffPt = new TCanvas("Pho_reco_eff_pt","RecoEffPt_compare", 800, 600);
  PhoRecoEffPt->cd();
  reco_eff_pt_8TeV.DrawCopy();
  reco_eff_pt_14TeV.DrawCopy("same");
  TLegend* L2 = new TLegend(x1,y1,x2,y2);
  L2->AddEntry(&reco_eff_pt_8TeV,  "8TeV",    "L");
  L2->AddEntry(&reco_eff_pt_14TeV, "14TeV",   "L");
  L2->Draw();
  PhoRecoEffPt->Update();
  PhoRecoEffPt->Write();



  outputFile->Close();


}
